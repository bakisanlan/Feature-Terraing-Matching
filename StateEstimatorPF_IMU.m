classdef StateEstimatorPF_IMU < handle
    %StateEstimatorPF Manages a state estimator that is based on terrain
    %matching using standard Particle Filter (PF). This is class is written
    %on purpose on TAN(Terrain Aided Navigation), so for general problem
    % MATLAB built-in Particle Filter can be used.


    % The state estimator recevies point cloud measurements.
    %
    %   Each PF particle is aware of its world frame x-y position namely
    %   (xw,yw). This position represents a hypothetical x-y position of
    %   the aircraft whose ground track location is being estimated.
    %
    %   Given a point cloud measurement (Xm,Ym,Zm) represented in the
    %   sensor's frame, the likelihood of an associated location (xw,yw) is
    %   evaluated using one of two approaches:
    %       - Each particle raycasts from its current reference map
    %         position (xw,yw) to get point cloud (Xp,Yp,Zp) represented in
    %         the particle's frame. Elevations Zp are then compared with
    %         Zm. Currently, the (Xp,Yp) information is not compared with
    %         (Xm,Ym).
    %       - Each particle uses the (Xm,Ym) pattern and its current
    %         reference map position (xw,yw) to directly read Zp elevation
    %         values at (xw,yw)+(Xm,Ym) to get a point cloud (Xm,Ym,Zp).
    %         The Zp elevations are then compared with Zm. Thus no
    %         raycasting is performed.
    %
    %       Note that the Zp values using both methods are not necessarily
    %       equal numerically.
    %

    % At each filtering iteration, PF particles move through process model
    % then collect measurement through above methods, then particles are
    % weighted based on likelihood that is found by pdf value of
    % MAE(Mean Absolute Error) over on normal distribution of mean is 0
    % and std is standard deviation of measurement device.
    % likelihood  = L(MAE value | N(mean = 0, std = alt_std). Multimodal
    % distribution is collected, after assigned every particles weight.
    % Estimation result can be obtained by using weighted mean of particles.


    properties

        priorPose   % Best prior guess of the pose in world frame
        %   positionLiDAR     x, y, z        (ENU world frame)
        %   orientationLiDAR  psi-theta-phi  (Yaw angle of sensor frame w.r.t. to world frame)

        hReferenceMapScanner    % Scanner to be used by the estimator

        % Particle Filter Variables
        N                  %Number particles
        state_bounds       %Particles states bounds
        alt_std            %Altitude standart error of measurement device
        process_std        %Process error of kinematic model of particles
        particles          %Particles states information
        exp_rate           %Exploration rate of PF algorithm is bounded on [0,1]
        dt                 %Delta time of particles prediction
        weights            %Weight of each particle
        meann              %Estimation result of PF at current time step
        count_est          %Count of estimation on simulation
        var                %Estimation standard deviation of PF at current time step
        radar_Z            %Store all radar elevation measurement on simulation in world frame
        particles_pc       %Store particles measurement point cloud data at current time step
        MAE_particles      %Mean Absolute Error(MAE) value for each particle w.r.t. aircraft measurement
        MAE_particles_hist %Store all MAE values of each particle on simulation
        batch_size         %Number of measurement will interpret as information to MAE process.
        batch_n_Part       %Number of element of batch is collected
        MAE_Batch_Part     %MAE value at each element of batch for each particle
        old_index          %Indexes before last resample
        part_radar_pattern_world 
        imp_indexes
        n_random_particles
        terrain_uniqueness
        part_max_w
    end
    %% PF Functions
    methods (Access = public)
        function self = StateEstimatorPF_IMU(N,state_bounds,alt_std,dt,varargin)
            % Defining properties
            if nargin < 4
                error('Not enough input arguments are given');

            elseif nargin == 4
                self.exp_rate = 0;
                self.batch_size = 1;
            elseif nargin == 5
                self.exp_rate = varargin{1};
                self.batch_size = 1;
            elseif nargin == 6
                self.exp_rate = varargin{1};
                self.batch_size = varargin{2};
            else
                error('Too much input arguments are given')
            end

            self.N = N;
            self.state_bounds = state_bounds;
   
            self.alt_std = alt_std;
            self.dt = dt;
            self.batch_n_Part = zeros(self.N,1);
            self.process_std = [10 10 1 0 0 0];
            % self.process_std = [0.5 0.5 0.1 0.01 0.01 0.01];
            % self.process_std = [0.5 0.5 0.1 0 0 0];
            % self.process_std = zeros(1,6);
            self.old_index = 1:self.N;
            self.count_est = 1;
            self.n_random_particles = 0;

            % Create uniformly distributed particles around aircraft with
            % the span range
            self.particles = self.create_particles(self.N,state_bounds);
            self.weights = ones(self.N,1)/self.N;

        end

        function param_estimate = getEstimate(self,u,ground_truth,ptCloudRadar,flagRAYCAST, modelF)

            %getEstimate Returns the location estimate in world frame.
            %
            %   Takes as input partially known current Pose info, namely
            %   altitude and orientation above mean sea level Z, the
            %   control input, the measurement point cloud, and a method
            %   flag. Only the current X and Y location info is unknown and
            %   to be estimated.

            %Move to particles based on kinematic model with same aircraft
            %input
            self.particle_step(u,modelF);
            % self.particles(1,1:3) = ground_truth;

            %self.particles(:,7:9) = rot;

            %Finding MAE(mean absolute error) correlation of each
            %particle w.r.t. radar Z value
            self.find_MAE_particles(ptCloudRadar,flagRAYCAST);  % Nx81

            % update weight based on MAE
            self.update_weights()

            % get mean and var value and store
            self.estimate();
            param_estimate = [self.meann, self.var];

            % counting estimation times
            self.count_est = self.count_est + 1;
        end

    end

    %% Functions that is used in Class as subfunctions
    methods (Access = private)

        function particles = create_particles(self,n_partc,varargin)

            if nargin == 3
                state_bound = varargin{1};  %nx2
                n_state = length(state_bound);
                particles = unifrnd(repmat(state_bound(:,1)',n_partc,1), repmat(state_bound(:,2)',n_partc,1),n_partc,n_state);

            elseif nargin == 4
                mu = varargin{1}; %nx1
                std = varargin{2}; 
                n_state = length(mu);
                particles = normrnd(repmat(mu',n_partc,1), repmat(std',n_partc,1), n_partc, n_state);
            else
                error('Incorrect number of inputs')
            end
            % create particles uniform distribution near to the guess of aircraft
            % gps lost position

            %particles = zeros(n_partc, n_state); %3 vel , 3 pos, 3 orient
            % %uniform distribution


            % % %gaussian dstrubition
            % particles(:, 1) = normrnd(init_pos(1), span_partc(1)*0.25, [n_partc 1]);
            % particles(:, 2) = normrnd(init_pos(2), span_partc(1)*0.25, [n_partc 1]);

        end

        function particle_step(self, u, modelF)
            %Move the aircraft based on the provided input.
            %
            %   The kinematics are represented by
            %    dx/dt      = v * cos(psi);
            %    dy/dt      = v * sin(psi);
            %    dz/dt      = 0;
            %    dpsi/dt    = omega;
            %
            %    u = [v; omega]
            %rng(5,'twister')
            u = reshape(u, numel(u),1); % column vector
            noise = (randn(self.N, length(u)) .* self.process_std);
            u = u' + noise; % 200x6

            % Updated based on input

            % There are 2 different kinematical model. Heading is state
            % variable in the first model, heading is input variable
            % in the second model.
            if modelF == 1
                % dPose1_dt = u(:,1) .* cos(self.particles(:,4));
                % dPose2_dt = u(:,1) .* sin(self.particles(:,4));
                % dPose3_dt = 0;
                % dPose4_dt = u(:,2);
                % 
                % self.particles(:, 1) = self.particles(:, 1) + self.dt .* dPose1_dt;% + (randn(self.N, 1) * self.process_std);
                % self.particles(:, 2) = self.particles(:, 2) + self.dt .* dPose2_dt;% + (randn(self.N, 1) * self.process_std);
                % self.particles(:, 3) = self.particles(:, 3) + self.dt .* dPose3_dt;
                % self.particles(:, 4) = self.particles(:, 4) + self.dt .* dPose4_dt;

                %u(:,3) = zeros(self.N,1);
                % dPose_dt = u;
                % self.particles(:,1:3)   = dPose_dt(:,1:3);
                % self.particles(:,4:end) = self.particles(:,4:end) + self.dt .* dPose_dt;
                % 
                % self.particles(:,6) = 600*ones(self.N,1);
                
                % acc model
                % psi = self.particles(1,7); theta = self.particles(1,8); phi = self.particles(1,9);
                % bodyROTned = [
                %     cos(theta)*cos(psi)                             , cos(theta)*sin(psi)                             ,  -sin(theta);
                %     sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi), sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi), sin(phi)*cos(theta);
                %     cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi), cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi), cos(phi)*cos(theta)
                %                 ]';
                % 
                % trans_acc = u(:,1:3);
                % trans_acc_cell = num2cell(trans_acc,2);
                % 
                % rotat_vel = u(:,4:6);
                % rotat_vel_cell = num2cell(rotat_vel,2);
                % 
                % rotated_trans_acc = cellfun(@(x) bodyROTned*x',trans_acc_cell,'UniformOutput',false);
                % rotated_trans_acc = cell2mat(rotated_trans_acc);
                % rotated_trans_acc = reshape(rotated_trans_acc,3,[])';
                % g = [0 0 -9.809800000000001];
                % rotated_trans_acc = g - rotated_trans_acc;
                % 
                % rotated_rotat_vel = cellfun(@(x) bodyROTned*x',rotat_vel_cell,'UniformOutput',false);
                % rotated_rotat_vel = cell2mat(rotated_rotat_vel);
                % rotated_rotat_vel = reshape(rotated_rotat_vel,3,[])';
                % 
                % 
                % self.particles(:,4:6) = self.particles(:,4:6) + self.dt .* self.particles(:,1:3); %pos
                % self.particles(:,1:3) = self.particles(:,1:3) + self.dt .* rotated_trans_acc;     %vel
                % self.particles(:,7:9) = self.particles(:,7:9) + self.dt .* rotated_rotat_vel;     %rot

                %self.particles(:,6) = 600*ones(self.N,1);

                % position kinematic model

                self.particles(:,1:3) = self.particles(:,1:3) + self.dt .* u(:,1:3);
                self.particles(:,4:6) = u(:,4:6);
                

            else
                dPose1_dt = u(:,1) .* cos(u(:,2));
                dPose2_dt = u(:,1) .* sin(u(:,2));
                dPose3_dt = 0;

                self.particles(:, 1) = self.particles(:, 1) + self.dt .* dPose1_dt;% + (randn(self.N, 1) * self.process_std);
                self.particles(:, 2) = self.particles(:, 2) + self.dt .* dPose2_dt;% + (randn(self.N, 1) * self.process_std);
                self.particles(:, 3) = self.particles(:, 3) + self.dt .* dPose3_dt;
                self.particles(:, 4) = u(:,2);
            end
        end

        function update_weights(self)
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
            self.MAE_particles = cellfun(@mean,self.MAE_Batch_Part);
            self.MAE_particles_hist(:,self.count_est) = self.MAE_particles;

            % calcualte likelihood, then update prior as posterior.
            likelihood = (1/(self.alt_std*2.506628274631)) .* exp(-0.5 .* (self.MAE_particles./self.alt_std).^2);
            self.weights = self.weights .* likelihood;
            %self.weights = self.weights .* (1/(self.alt_std*sqrt(2*pi))) .* exp(-0.5 .* (self.MAE_particles./self.alt_std).^2);

            %normalizing weights after each update
            self.weights = self.weights + 1e-300; % preventing 0 weights value
            self.weights = self.weights ./ sum(self.weights);

            % % sorting particles based on weight on descending order
            % [~,I] = sort(self.weights,'descend');
            % self.weights = self.weights(I);
            % self.particles = self.particles(I,:);
            % self.MAE_particles = self.MAE_particles(I);

            % Below code for resampling criteria. There are many methods to
            % choose when we should resample but we choose N_eff < N/2 criteria
            % If condition is provided, we are resampling our particles
            % with Systematic resampling methods.
            if self.neff < self.N/2
            % if false
                indexes = self.resample_Systematic;
                self.resample(indexes)
            end
        end

        function find_MAE_particles(self,ptCloudRadar,flagRAYCAST)
            %This function will find Mean Absolute Error(MAE) value through
            %comparing of each particle measurement with aircraft measurement
            %by either using sliding method or raycasting method.

            % Updating reference map scanner via prior estimate as
            % asssuming altitude and attitudes are known.

            % particless = self.particles(:,4:end);
            % particless(:,4:6) = rot;

            [Xs,Ys,Zs, xw,yw,Zw_all]= self.hReferenceMapScanner.lookupElevation( ...
                                 self.particles,ptCloudRadar,flagRAYCAST);

            % self.terrain_uniqueness(self.count_est) = std(mean(Zw_all(:,6),2),0,1);
            % disp(self.terrain_uniqueness(self.count_est))

            for i=1:self.N
                % Xr_s = ptCloudRadar.s.Location(:,1);
                % Yr_s = ptCloudRadar.s.Location(:,2);
                % Zr_s = ptCloudRadar.s.Location(:,3);
                % 
                % Xss = Xs(i,:)';
                % Yss = Ys(i,:)';
                % Zss = Zs(i,:)';

                Zw = self.particles(i,3) - Zw_all(i,:)';
                Zr_w = 600 - ptCloudRadar.w.Location(:,3);
                % Zw =  Zw_all(i,:)';
                % Zr_w = ptCloudRadar.w.Location(:,3);

               
                % Decide the best way to treat NANs.
                idx = or(isnan(Zw),isnan(Zr_w));
                % Deleting NaNs
                Zw(idx) = [];
                Zr_w(idx) = [];

                % Xss(idx) = [];
                % Yss(idx) = [];
                % Zss(idx) = [];
                % 
                % Xr_s(idx) = [];
                % Yr_s(idx) = [];
                % Zr_s(idx) = [];

                %particle_pc_w(idx,:) = [];

                % TERCOM like historical correlation if batch_size ~=1
                % If batch_size ==1, below if condition will be
                % unnecessary but it will work.
                if (self.batch_n_Part(i) ~= self.batch_size)
                    % increase batch_n of ith particle
                    self.batch_n_Part(i) = self.batch_n_Part(i) + 1;
                    if ~isempty(Zr_w)
                        self.MAE_Batch_Part{i,1}(self.batch_n_Part(i)) = mean(abs(Zw - Zr_w),1);
                        %self.MAE_Batch_Part{i,1}(self.batch_n_Part(i)) = mean([mean(abs(Xss - Xr_s),1),mean(abs(Yss - Yr_s),1),mean(abs(Zss - Zr_s),1)]);

                    else
                        self.MAE_Batch_Part{i,1}(self.batch_n_Part(i)) = 999;
                    end
                else
                    % drop first column of batch in every iteration
                    self.MAE_Batch_Part{i,1}(1) = [];
                    if ~isempty(Zr_w)
                        % add new MAE to end of batch array
                        % self.MAE_Batch_Part{i,1}(self.batch_size) =  mean(abs(Zw - Zr_w),1);
                        self.MAE_Batch_Part{i,1}(self.batch_size) =  mean((Zw - Zr_w).^2,1);
                        %self.MAE_Batch_Part{i,1}(self.batch_n_Part(i)) = sqrt(mean(abs(Xss - Xr_s),1)^2 + mean(abs(Yss - Yr_s),1)^2 + mean(abs(Zss - Zr_s),1)^2);

                    else
                        self.MAE_Batch_Part{i,1}(self.batch_size) = 999;
                    end
                end

                % Storing particles point cloud measurement data in world
                % frame for visualization
                %self.particles_pc{i} = particle_pc_w;

                % Dealing empyt value with replacing NaN
                % if isempty(self.particles_pc{i})
                %     self.particles_pc{i} = NaN(1,3);
                % end

                % Storing Radar elevation value for every estimation history
                self.radar_Z{1,self.count_est} = Zr_w;
                self.part_radar_pattern_world{i,self.count_est} = {xw(i,:),yw(i,:),Zw_all(i,:)};


            end
        end

        function estimate(self)
            % Taking weighted mean and var of particles for estimation
            % returns mean and variance of the weighted particles
            pos = self.particles(1:self.N-self.n_random_particles,1:3);
            weights = self.weights(1:self.N-self.n_random_particles) / sum(self.weights(1:self.N-self.n_random_particles));
            self.meann = sum(pos .* weights,1) / sum(weights);
            self.var  = sum((pos - self.meann).^2 .* weights,1) / sum(weights);
            self.part_max_w = self.particles(1,:);
            %self.var = sqrt(var_vec(1)^2 + var_vec(2)^2);
        end

        function [mu , var] = estimate_solid(self)
            % Taking weighted mean and var of particles for estimation
            % returns mean and variance of the weighted particles
            pos = self.particles(1:self.N-self.n_random_particles,1:3);
            weights = self.weights(1:self.N-self.n_random_particles) / sum(self.weights(1:self.N-self.n_random_particles));
            mu = sum(pos .* weights,1) / sum(weights);
            var  = sum((pos - self.meann).^2,1) / self.N;
            %self.var = sqrt(var_vec(1)^2 + var_vec(2)^2);
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
            self.n_random_particles = 0;

            % Exploration phase
            % When resampling is activated, new samples are chosed by indexed
            % array that is produced by resampling method with some
            % 1-exploration rate. And remaining particles are created
            % uniform randomly in range of last mean estimated value with exploration rate.
            if (self.exp_rate ~= 0) && ~isempty(self.meann)
                n_random_particles = min(round(self.N*self.exp_rate),self.N - 1);
                % indexes = randsample(indexes, self.N - n_random_particles);
                indexes = indexes(1:(self.N - n_random_particles));
                % self.estimate()
                [mu , var] = self.estimate_solid;
                std  = sqrt(var);
                self.n_random_particles = n_random_particles;
                ref_var = 400;
                ref_std = 240;
                max_val = 300;
                % bounds_random_particles = [mean(1) - min((ref_var/var(1))*ref_std,max_val)       mean(1) + min((ref_var/var(1))*ref_std,max_val);
                %                            mean(2) - min((ref_var/var(2))*ref_std,max_val)       mean(2) + min((ref_var/var(2))*ref_std,max_val);
                %                            mean(3)   mean(3) ;
                %                            self.particles(1,4)  self.particles(1,4);
                %                            self.particles(1,5)  self.particles(1,5);
                %                            self.particles(1,6)  self.particles(1,6)];

                
                % rand_particles = self.create_particles(self.n_random_particles,bounds_random_particles);

                mu_part = [mu' ;
                          self.particles(1,4:6)'];

                std_part = [(ref_std ./ std(1:2)').*30 ; 
                            zeros(4,1)];

                rand_particles = self.create_particles(self.n_random_particles,mu_part,std_part);
                self.particles = [self.particles(indexes,:) ; rand_particles];
            else
                self.imp_indexes = indexes;
                self.particles = self.particles(indexes,:);
            end

            disp('---------------------------------------')
            disp('--------------Resampled----------------')
            disp('---------------------------------------')

            % reset batch and batch_n of resampled particle
            self.MAE_Batch_Part(:) = {[]};
            self.batch_n_Part(:) = 0;

            % Reset weights of particles
            self.weights = ones(length(self.particles),1)/length(self.particles);
        end

    end
end


