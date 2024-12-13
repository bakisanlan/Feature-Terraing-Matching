classdef StateEstimatorMPF < handle

    properties
        particles
        particles_cov
        weights
        likelihood
        N
        dt
        R
        n_nonlin
        n_lin
        Accelerometer
        Gyroscope
        KalmanFiltersState              %nxN
        KalmanFiltersCovariance         %nxnxN
        X
        P
        hReferenceMapScanner
        MAE_particles
        MAE_particles_hist
        count_est
        partc_radar_hist
        terrain_uniqueness
        partc_Zs_error
    end

    methods (Access = public)

        function self = StateEstimatorMPF(N,mu_part,std_part,mu_kalman,cov_kalman,dt,R)
            % rng(32,"twister");

            self.N  = N;
            self.dt = dt;
            self.R  = R;
            self.n_nonlin = length(mu_part);
            self.n_lin    = length(mu_kalman);
            self.count_est = 0;
            self.X = zeros(15,1);
            self.P = diag(zeros(15,1));


            % Create uniformly distributed particles around aircraft with
            % the span range
            self.MPF_initilize(mu_part,std_part,mu_kalman,cov_kalman);
        end

        function param_estimate = getEstimate(self,u,X_nom,UAVImage)
            % rng(32,"twister");

            self.count_est = self.count_est + 1;

            rotM_nom = quat2rotm(X_nom(7:10)');

            self.predict(u, rotM_nom);
            
            self.find_likelihood_particles(UAVImage,X_nom)
           
            self.update()

            param_estimate.State      = self.X;
            param_estimate.Covariance = self.P;
            
        end
    end


    methods (Access = private)

        function MPF_initilize(self,mu_part,std_part,mu_kalman,cov_kalman)

            mu = mu_part; %nx1
            std = std_part;
            self.particles = normrnd(repmat(mu,1,self.N), repmat(std,1,self.N), self.n_nonlin ,self.N);
            self.particles_cov = diag(std.^2);

            %initialize weights equally
            self.weights = ones(1,self.N)/self.N;

            %create kalman filters of particles
            self.KalmanFiltersState      = repmat(mu_kalman,  [1, self.N]);
            self.KalmanFiltersCovariance = repmat(cov_kalman, [1, 1, self.N]);
        end

        function predict(self, u, rotM_nom)

            n = self.n_nonlin;
            l = self.n_lin;

            acc_m  = -u(1:3);  %matlab accel reading should be negate
            gyro_m = u(4:6);

            acc_nom  = acc_m  + self.Accelerometer.ConstantBias'; % bias sign is plus due to negatiation of matlab
            % acc_nom(3) = 0;
            gyro_nom = gyro_m - self.Gyroscope.ConstantBias';

            acc_nom_skew = [ 0         -acc_nom(3)  acc_nom(2);
                            acc_nom(3)  0          -acc_nom(1);
                           -acc_nom(2)  acc_nom(1)  0];

            rotM_drot = self.exp_rot(gyro_nom*self.dt); %eq 78

            
            f_nonlin = eye(n);
            A_nonlin = [eye(n)*self.dt zeros(n,l-n)];

            f_lin    = zeros(l,n);
            A_lin    = [eye(n)   -rotM_nom*acc_nom_skew*self.dt  -rotM_nom*self.dt  zeros(n);
                        zeros(n)  rotM_drot'*self.dt              zeros(3)         -eye(n)*self.dt;
                        zeros(n)  zeros(n)                        eye(n)            zeros(n);
                        zeros(n)  zeros(n)                        zeros(n)          eye(3)         ];

            C = zeros(1,l);

            % Qn  = diag([0 0 0]);
            Qn  = diag([1 1 0]);
            Qnl = zeros(n,l);
            Ql  = diag([self.Accelerometer.BiasInstability.^2 * self.dt^2, ...
                        self.Gyroscope.BiasInstability.^2     * self.dt^2, ...
                        self.Accelerometer.NoiseDensity.^2    * self.dt, ...
                        self.Gyroscope.NoiseDensity.^2        * self.dt]);
            Q = [Qn   Qnl
                 Qnl' Ql];

            % Filters Propagation
            c = 3;
            for i = 1:self.N

                % Particle Filter Propagation
                particle_current = self.particles(:,i);
                % self.particles(:,i) = f_nonlin * self.particles(:,i) + A_nonlin * self.KalmanFiltersState(:,i) + (A_nonlin * self.particles_cov * A_nonlin).^1/c;
                self.particles(:,i) = f_nonlin * self.particles(:,i) + ...
                                      A_nonlin * self.KalmanFiltersState(:,i) + ...
                                      ((A_nonlin * self.KalmanFiltersCovariance(:,:,i) * A_nonlin' + Qn).^1/c) * randn(3,1);
                self.particles(3,i) = 0;
                
                % Kalman Filter Propagation
                A_lin_t = A_lin;% - Qnl'/Qn*A_nonlin;
                Ql_t    = Ql;   % - Qnl'/Qn*Qnl;

                N = A_nonlin*self.KalmanFiltersCovariance(:,:,i)*A_nonlin' + Qn; %non linear state covariance related
                L = A_lin_t*self.KalmanFiltersCovariance(:,:,i)*A_nonlin'/N;     %

                % Z = A_nonlin*self.KalmanFiltersState(:,i) + (A_nonlin * self.particles_cov * A_nonlin).^1/c;
                Z = self.particles(:,i) - f_nonlin * particle_current;

                % self.KalmanFiltersState(:,i)        = A_lin_t*self.KalmanFiltersState(:,i) + Qnl'/Qn*Z + L*(Z - A_nonlin*self.KalmanFiltersState(:,i));
                self.KalmanFiltersState(:,i)        = A_lin_t*self.KalmanFiltersState(:,i) + L*(Z - A_nonlin*self.KalmanFiltersState(:,i)) + f_lin*particle_current;
                self.KalmanFiltersState(4:6,i)      = zeros(3,1); % rotation is known
                self.KalmanFiltersCovariance(:,:,i) = A_lin_t*self.KalmanFiltersCovariance(:,:,i)*A_lin_t' + Ql_t - L*N*L';
            end
        end


        function find_likelihood_particles(self,UAVImage,X_nom)
            %This function will find Mean Absolute Error(MAE) value through
            %comparing of each particle measurement with aircraft measurement
            %by either using sliding method or raycasting method.

            % Updating reference map scanner via prior estimate as
            % asssuming altitude and attitudes are known.

            % particless = self.particles(:,4:end);
            % particless(:,4:6) = rot;

            % X_nom(3) = true_
            % h_INS = X_nom(3);
            pos_hypo = X_nom(1:3) + self.particles;
            % rot_hypo = quat2eul(X_nom(7:10)')' + self.KalmanFiltersState(7:9,:) + radar_tilt;
            % rot_hypo = repmat(true_rot' + radar_tilt,1,self.N);
            
            X_hypo = pos_hypo';
            % X_hypo(:,3) = true_pos(3);

            rot_hypo = quat2eul(X_nom(7:10)')' + self.KalmanFiltersState(7:9,:);
            yaw = rad2deg(rot_hypo(1));

            self.likelihood = self.hReferenceMapScanner.find_likelihood(X_hypo(:,1:2),yaw,UAVImage);
        end


        function update(self)
            % Update weight based on pdf value

            % Below code block is the core of Particle Filter algorithm.
            % For each particle, likelihood  that is found by pdf value of
            % MAE(Mean Absolute Error) over on normal distribution of mean is 0
            % and std is standard deviation of measurement device.
            % likelihood  = L(MAE value | N(mean = 0, std = alt_std). Multimodal
            %
            % After found likelihood value, we are updating our prior
            % estimate with likelihood value as posterior estimate(Bayes theorem)

            % find means of batch elements for finding final MAE
            %self.MAE_particles = cellfun(@mean,self.MAE_Batch_Part);
            % self.MAE_particles_hist(:,self.count_est) = self.MAE_particles;

            % calcualte likelihood, then update prior as posterior.
            % likelihood = (1/(self.R*2.506628274631)) .* exp(-0.5 .* (self.MAE_particles./self.alt_std).^2);
            % likelihood = exp()
            % E = self.MAE_particles - 
            % PF Mesaurement Update
            self.weights = self.weights .* self.likelihood';
            % self.weights = self.weights .* (1/(sqrt(self.R)*sqrt(2*pi))) .* exp(-0.5 .* (self.MAE_particles./sqrt(self.R)).^2);

            %normalizing weights after each update
            self.weights = self.weights + 1e-300; % preventing 0 weights value
            self.weights = self.weights ./ sum(self.weights);

            %Compute estimate for nonlinear states
            xn = self.particles;
            xn_est = sum(xn .* self.weights,2) / sum(self.weights); % dimension of xn
            Pn_est = sum((xn - xn_est).^2 .* self.weights,2) / sum(self.weights);

            % self.X(1:self.n_nonlin) = xn_est;
            % self.P(1:self.n_nonlin, 1:self.n_nonlin) = diag(Pn_est);
           
            % Below code for resampling criteria. There are many methods to
            % choose when we should resample but we choose N_eff < N/2 criteria
            % If condition is provided, we are resampling our particles
            % with Systematic resampling methods.
            if self.neff < self.N/2
            % if false
                indexes = self.resample_Systematic;
                %self.resample(indexes)
                self.particles               = self.particles(:,indexes);
                self.KalmanFiltersState      = self.KalmanFiltersState(:,indexes);
                self.KalmanFiltersCovariance = self.KalmanFiltersCovariance(:,:,indexes);
    
                % Reset weights of particles
                self.weights = ones(1,self.N)/self.N;


                disp('---------------------------------------')
                disp('--------------Resampled----------------')
                disp('---------------------------------------')

            end

            %Kalman Filters Measurement Update
            %Due to C is zeros no need to Update KF
            % for i = 1:self.N
            %     M = C*self.KalmanFiltersCovariance(:,i)*C' + R;
            %     K = self.KalmanFiltersCovariance(:,i)*C'/M;
            % 
            %     self.KalmanFiltersCovariance(:,i) = self.KalmanFiltersCovariance(:,i) - K*M*K';
            %     self.KalmanFiltersState(:,i) = self.KalmanFiltersState(:,i) + K*self.MAE_particles(i);
            % end
            xl = self.KalmanFiltersState;
            xl_est = sum(xl .* self.weights,2) / sum(self.weights); % dimension of xl, mean can be use
            aux_cov = (xl - xl_est) .* (xl - xl_est);
            aux_cov = bsxfun(@times, reshape(aux_cov, self.n_lin, 1, self.N), eye(self.n_lin));
            aux_w   = bsxfun(@times, reshape(self.weights, 1, 1, self.N), ones(self.n_lin));
            Pl_est = sum((self.KalmanFiltersCovariance + aux_cov) .* aux_w ,[3]) / sum(self.weights);

            % self.X(self.n_nonlin+1:self.n_lin+self.n_nonlin) = xl_est;
            % self.P(self.n_nonlin+1:self.n_lin+self.n_nonlin, self.n_nonlin+1:self.n_lin+self.n_nonlin) = diag(Pl_est);

            self.X = [xn_est ; xl_est];
            self.P = [diag(Pn_est) zeros(self.n_nonlin,self.n_lin) ;
                      zeros(self.n_lin,self.n_nonlin)  Pl_est];
        end

        function estimate(self)
            xn = self.particles;
            xn_est = sum(xn .* self.weights,2) / sum(self.weights); % dimension of xn
            Pn_est = sum((xn - xn_est).^2 .* self.weights,2) / sum(self.weights);

            xl = self.KalmanFiltersState;
            xl_est = sum(xl .* self.weights,2) / sum(self.weights);
            Pl_est = sum((xl - xl_est).^2 .* self.weights,2) / sum(self.weights);

            self.X = [xn ; xl];
            self.P = [Pn_est ; Pl_est];
        end

        function indx = resample_Systematic(self)
            % There are many resampling methods but we choose Systematic
            % resample that resample particles based on cumulative sum of
            % weights
            % [~,I] = sort(self.weights);
            % self.weights = self.weights(I);
            % self.particles = self.particles(I,:);
            % self.MAE_particles = self.MAE_particles(I);

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

        function resample(self,indexes)
            % Resampling from indexes that is produces by sampling methods
            self.particles = self.particles(:,indexes);

            disp('---------------------------------------')
            disp('--------------Resampled----------------')
            disp('---------------------------------------')


            % Reset weights of particles
            self.weights = ones(length(self.particles),1)/length(self.particles);
        end


        function quat = exp_quat(d_rot)
            % This function create a quaternion based on exponential map to
            % small angular deviation

            d_rot = reshape(d_rot,[],1); %force d_rot to be column vector
            norm_rot = norm(d_rot);

            qw = cos(norm_rot/2);
            qr = (d_rot/norm_rot) * sin(norm_rot/2);

            quat = [qw ; qr];
        end

        function rotM = exp_rot(self,d_rot)

            % This function create a rotation matrix based on exponential map to
            % small angular deviation, 
            %eq 78 of Quaternion kinematics for the error-state Kalman Filter

            d_rot = reshape(d_rot,[],1); %force d_rot to be column vector
            norm_rot = norm(d_rot);

            u_vec = d_rot / norm_rot;
            u_skew = [ 0        -u_vec(3)  u_vec(2);
                       u_vec(3)  0        -u_vec(1);
                      -u_vec(2)  u_vec(1)  0];

            rotM = eye(3)*cos(norm_rot) + u_skew * sin(norm_rot) + u_vec*u_vec'*(1-cos(norm_rot));
        end
    end
end
