classdef INSbot < handle

    properties(SetAccess = public)

        dt                 %delta time
        g                  %gravitational vector m/s2
        NomState           %Nominal States
        TrueKinState       %True Kinematic States of UAV : [position, oriantation(euler)
        TrueKinQuat        %True Quaternion of Kinematic Model of UAV
        TrueState          %True INS States of UAV
        AccBias            %Acceleration ias term
        GyroBias           %Gyroscope bias term
        imuFs              %IMU Sample rate
        IMU                %MATLAB IMU object for dead reckoning motion
        IMUTruth           %MATLAB IMU object for ground truth motion
        prevTrueVelocity   %Previous time step True Velocity of Kinematic Model
        IMUreading         %IMU Readings that is output of IMU object
        IMUTruthreading    %No bias+no errors  IMU reading that is output of IMU object
        inertiaROTbody     %Inertia frame to body frame rotation matrix  

    end

    methods (Access = public)

        function obj = INSbot(initialState, dt,ReferenceFrame,IMUtype)

            if ReferenceFrame == 'NED'
                obj.g = [0 ; 0 ;  9.81];
            elseif ReferenceFrame == 'ENU'
                % obj.g = [0 ; 0 ; -9.80665];
                obj.g = [0 ; 0 ; -9.81];
            end

            % INS Initialization
            obj.NomState = reshape(initialState, numel(initialState),1);
            obj.dt = dt;

            obj.TrueKinState = initialState(1:3);
            obj.TrueKinQuat  = initialState(7:10);

            % IMU initialization
            obj.imuFs = 1/obj.dt;
            obj.IMU      = imuSensor('accel-gyro',ReferenceFrame=ReferenceFrame,SampleRate=obj.imuFs);
            obj.IMUTruth = imuSensor('accel-gyro',ReferenceFrame=ReferenceFrame,SampleRate=obj.imuFs);

            % There are different types of IMU
            % IMU1 : Tactical Grade IMU
            % IMU2 : Commercial Grade IMU
            % IMU3 : Perfect IMU(there is no bias or noise)
            if IMUtype == 1 %HG4930CA51 IMU(Expensive)
                %'MeasurementRange',20*9.81, ...
                DataSheetAccel = accelparams( ...
                    'ConstantBias', 1.7*9.81*1e-4, ...
                    'BiasInstability', 0.025*9.81*1e-3, ...
                    'NoiseDensity', 0.03/60);

                %'MeasurementRange',200*(pi/180), ...
                DataSheetGyro = gyroparams( ...
                    'ConstantBias', 7*(pi/180)/3600, ...
                    'BiasInstability', 0.25*(pi/180)/3600, ...
                    'NoiseDensity', 0.04*(pi/180)/sqrt(3600));

            elseif IMUtype == 2 % SCH16T-K01(Cheap)
                %'MeasurementRange',8.15*9.81, ...
                DataSheetAccel = accelparams( ...
                    'ConstantBias', 4.75*9.81*1e-4, ...
                    'BiasInstability', 10*0.02*9.81*1e-3, ...
                    'NoiseDensity', 10*5*0.03/60);

                %'MeasurementRange',62.5*(pi/180), ...
                DataSheetGyro = gyroparams( ...
                    'ConstantBias', 240*(pi/180)/3600, ... %??
                    'BiasInstability', 10*0.5*(pi/180)/3600, ...
                    'NoiseDensity', 10*0.04*(pi/180)/sqrt(3600));

            else %Matlab Default IMU which is perfect IMU
                DataSheetAccel = accelparams();
                DataSheetGyro = gyroparams();
            end

            obj.IMU.Accelerometer = DataSheetAccel;
            obj.IMU.Gyroscope     = DataSheetGyro;

            obj.IMUTruth.Accelerometer      = accelparams();
            obj.IMUTruth.Gyroscope          = gyroparams();

            
        end

        function move(obj,u)

            %MOVE Moves the UAV kinematic model based on the provided input.
            %Inputs are u(1) = Velocity, u(2) = Bank angle
            %States are X = [pN, pE, pD, heading, flight path angle, bank angle] 

            g = norm(obj.g);  %gravitational constant
            ng = 1/cos(u(2)); %ng = 1/cos(bank) to get zero pitch rate

            % Kinematic model for UAV motion
            dx = [u(1)*cos(obj.TrueKinState(4))*cos(obj.TrueKinState(5));
                  u(1)*sin(obj.TrueKinState(4))*cos(obj.TrueKinState(5));
                 -u(1)*sin(obj.TrueKinState(5));
                  (g*ng/u(1))*(sin(u(2))/cos(obj.TrueKinState(5)));       % yaw rate
                  (g/u(1))*(ng*cos(u(2)) - cos(obj.TrueKinState(5)));     % pitch rate
                  0];                                                     % roll rate

            % Inertia to body rotation matrix
            obj.inertiaROTbody = quat2rotm(obj.TrueKinQuat')'; %for ENU we have to take transpose

            % Euler rate to body rate matrix
            eulerRate2bodyRate = [1  0                       -sin(obj.TrueKinState(5)); 
                                  0  cos(obj.TrueKinState(6)) sin(obj.TrueKinState(6))*cos(obj.TrueKinState(5));
                                  0 -sin(obj.TrueKinState(6)) cos(obj.TrueKinState(6))*cos(obj.TrueKinState(5))];

            % store angular velocities
            angular_vel_inertia = dx([6 5 4]);
            angular_vel_body    = eulerRate2bodyRate * angular_vel_inertia;
            
            % Rotation dynamics of kinematic equation
            obj.TrueKinQuat =  (quatmultiply(obj.TrueKinQuat' , (obj.exp_quat(angular_vel_body .* obj.dt))'))';
            eul = quat2eul(obj.TrueKinQuat');
            obj.TrueKinState([4 5 6]) = eul;

            % Positional dynamics of kinematic equation
            obj.TrueKinState([1 2 3]) = obj.TrueKinState([1 2 3])  + dx([1 2 3]) * obj.dt;

            % post move calculation for getting velocity after step to find
            % translational acceleration with using previous velocity.
            trans_vel_post = [u(1)*cos(obj.TrueKinState(4))*cos(obj.TrueKinState(5));
                              u(1)*sin(obj.TrueKinState(4))*cos(obj.TrueKinState(5));
                             -u(1)*sin(obj.TrueKinState(5))];

            % update sensor reading
            obj.updateSensorReadings([trans_vel_post ; angular_vel_inertia]);
        end

        function moveINS(obj)

            %%Bias+Noise IMU Move(Represent real flight dead reckoning
            %%motion) which is Nominal State
            RotM = quat2rotm(obj.NomState(7:10)');
            acc_m  = -obj.IMUreading{1};     % measured acceleration, matlab accel reading should be negate
            % gyro_m =  obj.IMUreading{2};   % measured angular rate, use it when rotation contains errors
            gyro_m = obj.IMUTruthreading{2}; % measured angular rate, use it when rotation is known exactly

            acc_nom  = acc_m'  + obj.NomState(11:13);   %bias sign is plus due to negatiation of acc vector on MATLAB built in IMU sensor
            % gyro_nom = gyro_m' - obj.NomState(14:16); % use it when rotation contains errors
            gyro_nom = gyro_m' - 0;                     % use it when rotation is known exactly


            % Nominal state kinematics : [position, velocity, quaternion, acc, bias, gyro bias]
            V = obj.NomState(4:6);              
            V(3) = 0;   % use it when there is no any vertical velocity
            obj.NomState(1:3)   = obj.NomState(1:3) + obj.dt * (V +  0*0.5 * (RotM * acc_nom + obj.g) * obj.dt);
            obj.NomState(4:6)   = obj.NomState(4:6) + obj.dt * (RotM * acc_nom + obj.g);
            obj.NomState(7:10)  = quatmultiply(obj.NomState(7:10)' , (obj.exp_quat(gyro_nom * obj.dt))');
            obj.NomState(11:13) = obj.NomState(11:13);
            obj.NomState(14:16) = obj.NomState(14:16);

            %%Ground Truth Move(Represent ground truth UAV motion to track
            %%real states). There are no any noise or biased in that
            %%acceleration or angular rates
            RotM = quat2rotm(obj.TrueState(7:10)');
            acc_m  = -obj.IMUTruthreading{1};  % matlab accel reading should be negate
            gyro_m =  obj.IMUTruthreading{2};

            acc_true  = acc_m'  + 0; %bias sign is plus due to negatiation of acc vector on matlab
            gyro_true = gyro_m' - 0;

            % Ground Truth State Kinematics due to there is no error in 
            V = obj.TrueState(4:6);
            obj.TrueState(1:3)   = obj.TrueState(1:3) + obj.dt * (V +  0*0.5 * (RotM * acc_true + obj.g) * obj.dt);
            obj.TrueState(4:6)   = obj.TrueState(4:6) + obj.dt * (RotM * acc_true + obj.g);
            obj.TrueState(7:10)  = quatmultiply(obj.TrueState(7:10)' , (obj.exp_quat(gyro_true * obj.dt))');
            obj.TrueState(11:13) = obj.TrueState(11:13);
            obj.TrueState(14:16) = obj.TrueState(14:16);

        end

        function CorrectedState = correctINS(obj,dx)
            %This function corrects the shifted INS with MPF error states
            %solution 'dx'.
            %Correction is made by addition except quaternion is
            %multiplicalition for correction.
   
            % Multiplicative correction term on quaternion from error
            dq = obj.exp_quat(dx(7:9));
            qt = quatmultiply(obj.NomState(7:10)',dq')';  %corrected quaternion

            dx_noRot   = [dx(1:6)           ; dx(10:15)];            %error states without rotation information
            Xnom_noRot = [obj.NomState(1:6) ; obj.NomState(11:16)];  %nominal states without rotation information

            CorrectedState_noRot = Xnom_noRot + dx_noRot;            %correction of states without rotation
            CorrectedState       = [CorrectedState_noRot(1:6) ; qt ; CorrectedState_noRot(7:12)]; %getting corrected states with rotation
        end
    end

    methods (Access = private)

        function quat = exp_quat(obj,d_rot)
            % This function create a quaternion based on exponential map to
            % small angular deviation

            norm_rot = norm(d_rot);

            qw = cos(norm_rot/2);
            qr = (d_rot/norm_rot) * sin(norm_rot/2);
            if any(isnan(qr))
                qr = zeros(3,1);
            end

            quat = [qw ; qr];
        end

        function updateSensorReadings(obj, trueVelocity)
            %This function updates IMU readings with using translational
            %and angular velocity information from kinematic UAV dynamics

            % Get the true position, velocity, and acceleration in local 
            % navigation frame.
            trueTranslationalVelocity = trueVelocity(1:3);
            trueAngularVelocity       = trueVelocity(4:6); % the order of angular velocity should be X,Y,Z

            % Calculating acceleration value with finite differentiation
            % using previous and current velocity information
            trueAcceleration = (trueTranslationalVelocity - obj.prevTrueVelocity(1:3))/obj.dt;
            obj.prevTrueVelocity = trueVelocity;

            % Update MATLAB built in IMU model to get body reference value 
            [accelReading_body, gyroReading_body] = obj.IMU(trueAcceleration',trueAngularVelocity',obj.inertiaROTbody);
            obj.IMUreading  = {accelReading_body, gyroReading_body};

            % To track ground truth information, we track ground truth body
            % frame reading with using zero biased and zero noised IMU
            % model which is obj.IMUTruth
            [accelReading_body_truth, gyroReading_body_truth] = obj.IMUTruth(trueAcceleration',trueAngularVelocity',obj.inertiaROTbody);
            obj.IMUTruthreading  = {accelReading_body_truth, gyroReading_body_truth};

            % Call 'moveINS' function to step INS equations with body frame
            % inputs
            obj.moveINS

        end

    end

end