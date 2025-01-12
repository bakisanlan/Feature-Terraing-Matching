classdef StateEstimatorMPF < handle
    %%This is the class structure for Marginalized Particle Filter(MPF)
    %Main property of MPF, states of filter is divided two part: nonlinear
    %and linear. Nonlinear states represent nonlineer Monte Carlo filtering
    %prosedure. Mainly particles measure update based on measurement
    %function M(pos particles, rotation particles) that is highly
    %nonlinear. There is no any linear terms associate with measurement
    %update of particles. Because of rotation errors are negligible we can
    %assume that nonlinear states are solely positional error states. Other
    %states velocity, orientation, accelerometer bias and gyroscope bias
    %can be thought as linear and standard Kalman Filter update equation
    %can be used. Hence, with using MPF we combine both Monte Carlo
    %nonlinear filtering and linear Kalman Filtering. 

    properties
        N                         % Number of particles
        particles                 % (3 x N) particles position matrix        
        %particles_cov            % 
        weights                   % (1 x N) weights matrix of particles 
        likelihood                % (1 x N) likelihood matrix of particles
        dt                        % Delta time
        R                         % Measurement error for likelihood calculation
        n_nonlin                  % Number of nonlinear states which is 3
        n_lin                     % Number of linear states which is 12
        Accelerometer             % Accelerometer object of IMU   
        Gyroscope                 % Gyroscope object of IMU
        KalmanFiltersState        %(n_lin x N) Kalman Filter State matrix that store each particle Linear states
        KalmanFiltersCovariance   %(n_lin x n_lin x N) Kalman Filter State Covarianca matrix that store each particle Linear covariance states
        X                         % (15 x 1) MPF combined linear and nonlinear estimated state
        P                         % (15 x 15) MPF combined linear and nonlinear estimated state covariance
        hDataBaseScanner          % DataBaseScanner object that get measurement from offline database
        count_est                 % Count how much estimation made
    end

    methods (Access = public)

        function self = StateEstimatorMPF(N,mu_part,std_part,mu_kalman,cov_kalman,dt,R)

            % Set MPF variables 
            self.N  = N;
            self.dt = dt;
            self.R  = R;
            self.n_nonlin = length(mu_part);
            self.n_lin    = length(mu_kalman);
            self.count_est = 0;
            self.X = zeros(15,1);
            self.P = diag(zeros(15,1));


            % Create gaussian distributed particles around aircraft with
            % the meand and std values
            self.MPF_initilize(mu_part,std_part,mu_kalman,cov_kalman);
        end

        function param_estimate = getEstimate(self,u,X_nom,UAVImage)
            % This function is estimation loop of MMP that consist of
            % error states prediction, likelihood calculation of each
            % particles, measurement update of filter and give estimation of
            % states as output 'param_estimate'

            % Track estimation count
            self.count_est = self.count_est + 1;

            % Prediction step of MPF states
            self.predict(u, X_nom);
            
            % Find likelihood of particles through database
            self.find_likelihood_particles(UAVImage,X_nom)
           
            % Measurement update of MPF
            self.update()

            % Get estimated states and covariance matrix of MPF solution 
            self.estimate()
            
            % Resample particles based on weights
            self.resample()

            % Get estimation parameters
            param_estimate.State      = self.X;
            param_estimate.Covariance = self.P;
        end
    end


    methods (Access = private)

        function MPF_initilize(self,mu_part,std_part,mu_kalman,cov_kalman)

            % Inizilize particles with gaussian distribution [mu_part, std_part] 
            self.particles = normrnd(repmat(mu_part,1,self.N), repmat(std_part,1,self.N), self.n_nonlin ,self.N);

            % Initialize weights of particles equally as 1/N
            self.weights = ones(1,self.N)/self.N;

            % Create linear Kalman Filters of each particles
            self.KalmanFiltersState      = repmat(mu_kalman,  [1, self.N]);
            self.KalmanFiltersCovariance = repmat(cov_kalman, [1, 1, self.N]);
        end

        function predict(self, u, X_nom)
            %This function propagate error states kinematic equation given
            %IMU u vector that is [body acceleration(3 x 1) ; body angular rate(3 x 1)]

            n = self.n_nonlin;
            l = self.n_lin;

            % Get measured IMU outputs
            acc_m  = -u(1:3);  %matlab accel reading should be negate
            gyro_m =  u(4:6);

            % Get nominal IMU outputs with subtracting bias to measured ones
            acc_nom  = acc_m  + X_nom(11:13); % bias sign is plus due to negatiation of MATLAB IMU on acceleration
            % acc_nom(3) = 0;                 % Use it when there is no velocity change in Z direction
            gyro_nom = gyro_m - X_nom(14:16);


            % Skew symetric matrix of nominal acceleration 
            acc_nom_skew = [ 0         -acc_nom(3)  acc_nom(2);
                            acc_nom(3)  0          -acc_nom(1);
                           -acc_nom(2)  acc_nom(1)  0];

            % Rotation matrix that of nominal gyro rates cause
            rotM_drot = self.exp_rot(gyro_nom*self.dt); %eq 78

            % Nominal rotation matrix body to inertial
            rotM_nom = quat2rotm(X_nom(7:10)');

            %%MPF error state prediction equations are listed below. For
            %%the reference you can check :'Marginalized particle filters 
            %%for mixed linear/nonlinear state-space models' 
            %%https://ieeexplore.ieee.org/document/1453762

            % particles(t+1) = f_nonlin*particles(t) + A_nonlin*KalmanFilterState(t) + (A_nonlin*KalmanFilterCovariance(t)*A_nonlin' + Qn)^(1/c)*randn(3,1)
            % KalmanFilterState(t+1) = A_lin_t*KalmanFiltersState(t) + L*(Z- A_nonlin*KalmanFiltersState(t) + f_lin*particles(t)
            % KalmanFiltersCovariance(t+1) =
            % A_lin_t*KalmanFiltersCovariance(t)*A_lin_t' + Ql_t - L*N*L'   
            % N = A_nonlin*KalmanFiltersCovariance(t)*A_nonlin' + Qn %Dont confuse this N is not related with number of particles, this is a Kalman Optimality terms
            % L = A_lin_t*KalmanFiltersCovariance(t)*A_nonlin'*inv(N)
            % Z = particles(t+1) - f_nonlin * particles(t);
            % A_lin_t = A_lin - Qnl'/Qn*A_nonlin;
            % Ql_t    = Ql    - Qnl'/Qn*Qnl;

            % Nonlinear states related matrix
            f_nonlin = eye(n);
            A_nonlin = [eye(n)*self.dt zeros(n,l-n)];

            % Linear states related matrix
            f_lin    = zeros(l,n);
            A_lin    = [eye(n)   -rotM_nom*acc_nom_skew*self.dt  -rotM_nom*self.dt  zeros(n);
                        zeros(n)  rotM_drot'                      zeros(3)         -eye(n)*self.dt;
                        zeros(n)  zeros(n)                        eye(n)            zeros(n);
                        zeros(n)  zeros(n)                        zeros(n)          eye(3)         ];

            % Measurement function linear states related matrix 
            C = zeros(1,l); %there is no relation between linear states and measurement function

            % Nonlinear states process noise covariance matrix, which can 
            % be thought as particles process noise in positional states.
            % When Qn is high more exploration will be made by particles, 
            % when Qn is low particles can be missconverge somewhere.
            % Qn  = diag([0 0 0]); %most of case do not use zero process noise
            Qn  = diag([1 1 0]);
            
            % Covariance matrix between nonlinear and linear states which
            % is zero for every value
            Qnl = zeros(n,l);

            % Linear states covariance matrix which indicate IMU sensor
            % noise in linear states
            Ql  = diag([self.Accelerometer.NoiseDensity.^2               * self.dt, ...
                        self.Gyroscope.NoiseDensity.^2                   * self.dt, ...
                        ((2*self.Accelerometer.BiasInstability.^2)/3600) * self.dt, ...
                        ((2*self.Gyroscope.BiasInstability.^2)/3600)     * self.dt]); 

            % Augmented covariance matrix of full states
            Q = [Qn   Qnl
                 Qnl' Ql];

            % MPF Error States Propagation loop for each particles. Because
            % of every particles have a Kalman Filter inside, we loop
            % propagation equation through each particles.
            c = 2; % Cholesky Decomposition contant 
            for i = 1:self.N

                % Nonlinear States Propagation
                particle_current = self.particles(:,i); % particles in t time
                self.particles(:,i) = f_nonlin * self.particles(:,i) + ...
                                      A_nonlin * self.KalmanFiltersState(:,i) + ...
                                      ((A_nonlin * self.KalmanFiltersCovariance(:,:,i) * A_nonlin' + Qn).^1/c) * randn(3,1) * 0;
                self.particles(3,i) = 0; %use it when altitude is known exactly
                
                % Linear States Propagation
                A_lin_t = A_lin;% - Qnl'/Qn*A_nonlin;
                Ql_t    = Ql;   % - Qnl'/Qn*Qnl;
                N = A_nonlin*self.KalmanFiltersCovariance(:,:,i)*A_nonlin' + Qn; %non linear state covariance related
                L = A_lin_t*self.KalmanFiltersCovariance(:,:,i)*A_nonlin'/N;     
                Z = self.particles(:,i) - f_nonlin * particle_current;
                self.KalmanFiltersState(:,i)        = A_lin_t*self.KalmanFiltersState(:,i) + L*(Z - A_nonlin*self.KalmanFiltersState(:,i)) + f_lin*particle_current;
                self.KalmanFiltersState(4:6,i)      = zeros(3,1); % rotation is known
                self.KalmanFiltersCovariance(:,:,i) = A_lin_t*self.KalmanFiltersCovariance(:,:,i)*A_lin_t' + Ql_t - L*N*L';
            end
        end


        function find_likelihood_particles(self,UAVImage,X_nom)
            %This function will find likelihood value for each particles
            %with making image matching between UAVImage and each particles
            %image taken from database.

            % Hypotetical position of particles can be found by adding
            % nominal position of INS with position error that is
            % calculated by MPF.
            pos_hypo = X_nom(1:3) + self.particles;
           
            % Use below hypotetical rotation similar to hypotetical
            % position when rotatio is not known. Lets assume that we know
            % orientation of UAV exactly so there is no any error in
            % orientation. In this version we are trying to estimate only
            % position on 2D XY with using truth rotation information, also
            % we assume that there is just only yaw deflection. No bank or
            % pitching is allowed. In the next version we will estimate
            % both position and orientation.
            qdt  = self.exp_quat(self.KalmanFiltersState(4:6,:));
            qcorrected =  quatmultiply(X_nom(7:10)', qdt');
            rot_hypo = quat2eul(qcorrected)';
            yaw = rad2deg(rot_hypo(1)); 

            X_hypo = pos_hypo';
            %X_hypo = [pos_hypo ; rot_hypo]'; % Use it when rotation is not known.
            % X_hypo(:,3) = true_pos(3);      % Use it when altitude is known

            % Find likelihood of each particles using XY position of
            % particles(assuming altitude is known) and yaw rotation is
            % known and pitch and roll is zero.
            self.likelihood = self.hDataBaseScanner.find_likelihood(X_hypo(:,1:2),yaw,UAVImage);
        end


        function update(self)
            % This function block is the core of Particle Filter algorithm.
            % For each particle, posterior weight are updated with 
            % likelihood value that is found by amount of matched point 
            % between UAV image and particle.
            % After found likelihood value, we are updating our prior
            % estimate with likelihood value as posterior estimate(Bayes theorem)

            % PF Mesaurement Update
            self.weights = self.weights .* self.likelihood';

            %normalizing weights after each update
            self.weights = self.weights + 1e-300; % preventing 0 weights value
            self.weights = self.weights ./ sum(self.weights);
        end

        function estimate(self)
            %This function calculate estimated states and estimated
            %covariance matrix of MPF solution with using weights.
            xn = self.particles;
            xn_est = sum(xn .* self.weights,2) / sum(self.weights); % dimension same with xn
            Pn_est = sum((xn - xn_est).^2 .* self.weights,2) / sum(self.weights);


            %Kalman Filters Measurement Update
            %Due to C is zeros no need to Update KF
            % for i = 1:self.N
            %     M = C*self.KalmanFiltersCovariance(:,i)*C' + R;
            %     K = self.KalmanFiltersCovariance(:,i)*C'/M;
            % 
            %     self.KalmanFiltersCovariance(:,i) = self.KalmanFiltersCovariance(:,i) - K*M*K';
            %     self.KalmanFiltersState(:,i) = self.KalmanFiltersState(:,i) + K*self.MAE_particles(i);
            % end

            %Linear states mean and covariance compuatation
            xl = self.KalmanFiltersState;
            xl_est = sum(xl .* self.weights,2) / sum(self.weights); % dimension is same with xl
            %Diagonal covariance matrix calculation through xl_est and KalmanFiltersState
            aux_cov = (xl - xl_est) .* (xl - xl_est);
            aux_cov = bsxfun(@times, reshape(aux_cov, self.n_lin, 1, self.N), eye(self.n_lin)); %bsxfun function is an operator make element-wise @times operation A and B matrix
            aux_w   = bsxfun(@times, reshape(self.weights, 1, 1, self.N), ones(self.n_lin));
            Pl_est = sum((self.KalmanFiltersCovariance + aux_cov) .* aux_w ,[3]) / sum(self.weights);

            % Create combined estimated error state X array and estimated
            % covariance matrix of error state.
            self.X = [xn_est ; xl_est];
            self.P = [diag(Pn_est) zeros(self.n_nonlin,self.n_lin) ;
                      zeros(self.n_lin,self.n_nonlin)  Pl_est];
        end

        function indx = resample_Systematic(self)
            % There are many resampling methods but we choose Systematic
            % Resample that resample particles based on cumulative sum of
            % weights
            Q = cumsum(self.weights);
            indx = zeros(1);
            T = linspace(0,1-1/self.N,self.N) + rand(1)/self.N;
            T(self.N+1) = 1;
            i=1;
            j=1;

            while (i<=self.N)
                if (T(i)<Q(j))
                    indx(i)=j;
                    i=i+1;
                else
                    j=j+1;
                end
            end
        end

        function n_eff = neff(self)
            % Calculating N_eff for triggering of resampling
            n_eff = 1 ./ sum((self.weights).^2);
        end

        function resample(self)
            % This function for resampling of MPF. There are many methods to
            % choose. Resample should be activated when N_eff < N/2 criteria
            % If condition is provided, we are resampling our particles
            % with Systematic resampling methods.

            if self.neff < self.N/2

                indexes = self.resample_Systematic;

                % Resampling from indexes that is produces by sampling methods
                self.particles               = self.particles(:,indexes);
                self.KalmanFiltersState      = self.KalmanFiltersState(:,indexes);
                self.KalmanFiltersCovariance = self.KalmanFiltersCovariance(:,:,indexes);
    
                % Reset weights of particles to 1/N
                self.weights = ones(1,self.N)/self.N;

                % LINEAR ERROR STATE RESET
                % For now, we use open-loop manner, if we choose select closed
                % loop manner we have to use below reset command
                % self.KalmanFiltersState = repmat(zeros(12,1),  [1, self.N]);
    
                % Indicate resampling is made
                disp('---------------------------------------')
                disp('--------------Resampled----------------')
                disp('---------------------------------------')

            end
        end


        function quat = exp_quat(self,d_rot)
            % This function create a quaternion based on exponential map to
            % small angular deviation

            % d_rot is 3xM matrix

            norm_rot = vecnorm(d_rot,2,1);
            
            qw = cos(norm_rot/2);
            qr = (d_rot./norm_rot) .* sin(norm_rot/2);

            % dealt with NaN values which mean zero rotation
            idx_nan = ~isfinite(qr(1,:));
            qr(:,idx_nan) = zeros(3,sum(idx_nan));
            
            quat = [qw ; qr];

        end

        function rotM = exp_rot(self,d_rot)
            % This function create a rotation matrix based on exponential map to
            % small angular deviation, 
            %eq 78 of Quaternion kinematics for the error-state Kalman Filter

            d_rot = reshape(d_rot,[],1); %force d_rot to be column vector
            norm_rot = norm(d_rot);

            if norm_rot == 0

                rotM = eye(3);
                return;
            end

            u_vec = d_rot / norm_rot;
            u_skew = [ 0        -u_vec(3)  u_vec(2);
                       u_vec(3)  0        -u_vec(1);
                      -u_vec(2)  u_vec(1)  0];

            rotM = eye(3)*cos(norm_rot) + u_skew * sin(norm_rot) + u_vec*u_vec'*(1-cos(norm_rot));
        end
    end
end
