classdef INSbot < handle

    properties(SetAccess = public)

        dt
        g
        NomState
        TrueKinState
        TrueState
        TrueKinQuat
        AccBias
        GyroBias
        imuFs
        kinFs
        KINRateIMU
        IMU
        IMUTruth
        KINPerIMU
        nedROTbody
        prevTrueVelocity
        IMUreading
        IMUTruthreading
        insdt
        kindt
        %NomPrevState

    end

    methods (Access = public)

        function obj = INSbot(initialState, INSTs, KINTs,ReferenceFrame,IMUtype)

            if ReferenceFrame == 'NED'
                obj.g = [0 ; 0 ;  9.81];
            elseif ReferenceFrame == 'ENU'
                % obj.g = [0 ; 0 ; -9.80665];
                obj.g = [0 ; 0 ; -9.81];
            end

            % INS Initialization
            obj.NomState = reshape(initialState, numel(initialState),1);
            obj.insdt = INSTs;
            obj.kindt = KINTs;
            obj.TrueKinState = initialState(1:3);
            obj.TrueKinQuat  = initialState(7:10);

            % IMU initialization
            obj.imuFs = 1/INSTs;
            obj.kinFs = 1/KINTs;
            obj.KINRateIMU = obj.kinFs/obj.imuFs;
            obj.KINPerIMU = 0;

            obj.IMU      = imuSensor('accel-gyro',ReferenceFrame='ENU',SampleRate=obj.imuFs);
            obj.IMUTruth = imuSensor('accel-gyro',ReferenceFrame='ENU',SampleRate=obj.imuFs);

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
                    'BiasInstability', 100*0.02*9.81*1e-3, ...
                    'NoiseDensity', 5*0.03/60);

                %'MeasurementRange',62.5*(pi/180), ...
                DataSheetGyro = gyroparams( ...
                    'ConstantBias', 240*(pi/180)/3600, ... %??
                    'BiasInstability', 100*0.5*(pi/180)/3600, ...
                    'NoiseDensity', 5*0.04*(pi/180)/sqrt(3600));

            else %Matlab Default IMU
                DataSheetAccel = accelparams();
                DataSheetGyro = gyroparams();
            end

            obj.IMU.Accelerometer = DataSheetAccel;
            obj.IMU.Gyroscope     = DataSheetGyro;

            obj.IMUTruth.Accelerometer      = accelparams();
            obj.IMUTruth.Gyroscope          = gyroparams();

            
        end

        function move(obj,u,modelF)

            %MOVE Moves the aircraft based on the provided input.
            %
            %   The kinematics are represented by
            %    dx/dt      = v * cos(psi);
            %    dy/dt      = v * sin(psi);
            %    dz/dt      = 0;
            %    dpsi/dt    = omega;
            %
            %    u = [v; omega]

            % obj.nedROTbody = [cos(obj.KinState(4)) sin(obj.KinState(4)) 0
            %                  -sin(obj.KinState(4)) cos(obj.KinState(4)) 0
            %                   0                    0                    1];

            u = reshape(u, numel(u),1); % column vector

            % Updated based on input

            % There are 2 different kinematical model. Heading is state
            % variable in the first model, heading is input variable
            % in the second model.
            if modelF == 1
                % f =[0;
                %     0;
                %     0;
                %     0];
                %
                % g =[cos(obj.Pose(4)) 0;
                %     sin(obj.Pose(4)) 0;
                %     0                0;
                %     0                1];

                f =[0;
                    0;
                    0];

                g =[cos(obj.TrueKinState(4)) 0;
                    sin(obj.TrueKinState(4)) 0;
                    0                    0;];

                obj.nedROTbody = quat2rotm(obj.TrueKinQuat')';

                obj.TrueKinQuat =  (quatmultiply(obj.TrueKinQuat' , (obj.exp_quat([0 ; 0 ; u(2)] .* obj.kindt))'))';
                eul = quat2eul(obj.TrueKinQuat');


                 % quatmultiply(obj.NomState(7:10)' , (obj.exp_quat(gyro_nom * obj.insdt))');

                obj.TrueKinState(4) = eul(1);

                
                dPose_dt = f + g * (u);
                % obj.Pose = obj.Pose  + dPose_dt *  obj.dt;
                obj.TrueKinState(1:3) = obj.TrueKinState(1:3)  + dPose_dt * obj.kindt;

            else
                dPose_dt  = [u(1)*cos(u(2)) ; u(1)*sin(u(2)) ; 0];
                obj.TrueKinState(1:3) = obj.TrueKinState(1:3) + dPose_dt * obj.kindt;
                obj.TrueKinState(4) = u(2);
            end

            % Post move updates
            g =[cos(obj.TrueKinState(4)) 0;
                sin(obj.TrueKinState(4)) 0;
                0                0;];
            obj.updateSensorReadings([g ; 0 1]*u);

        end

        function moveINS(obj)

            %%IMU Move
            RotM = quat2rotm(obj.NomState(7:10)');
            acc_m  = -obj.IMUreading{1};  % matlab accel reading should be negate
            gyro_m =  obj.IMUreading{2};

            acc_nom  = acc_m'  + obj.NomState(11:13); %bias sign is plus due to negatiation of acc vector on matlab
            gyro_nom = gyro_m' - obj.NomState(14:16);

            V = obj.NomState(4:6);
            % V(3) = 0;

            % Nominal state kinematics
            obj.NomState(1:3)   = obj.NomState(1:3) + obj.insdt * (V +  0*0.5 * (RotM * acc_nom + obj.g) * obj.insdt);
            obj.NomState(4:6)   = obj.NomState(4:6) + obj.insdt * (RotM * acc_nom + obj.g);
            obj.NomState(7:10)  = quatmultiply(obj.NomState(7:10)' , (obj.exp_quat(gyro_nom * obj.insdt))');
            obj.NomState(11:13) = obj.NomState(11:13);
            obj.NomState(14:16) = obj.NomState(14:16);

            %%Ground Truth Move
            RotM = quat2rotm(obj.TrueState(7:10)');
            acc_m  = -obj.IMUTruthreading{1};  % matlab accel reading should be negate
            gyro_m =  obj.IMUTruthreading{2};

            acc_true  = acc_m'  + 0; %bias sign is plus due to negatiation of acc vector on matlab
            gyro_true = gyro_m' - 0;

            V = obj.TrueState(4:6);
            % V(3) = 0;

            % Nominal state kinematics
            obj.TrueState(1:3)   = obj.TrueState(1:3) + obj.insdt * (V +  0*0.5 * (RotM * acc_true + obj.g) * obj.insdt);
            obj.TrueState(4:6)   = obj.TrueState(4:6) + obj.insdt * (RotM * acc_true + obj.g);
            obj.TrueState(7:10)  = quatmultiply(obj.TrueState(7:10)' , (obj.exp_quat(gyro_true * obj.insdt))');
            obj.TrueState(11:13) = obj.TrueState(11:13);
            obj.TrueState(14:16) = obj.TrueState(14:16);
        end

        function CorrectedState = correctINS(obj,dx)
   
            % dq = obj.exp_quat(dx(7:9));
            dq = obj.exp_quat(dx([9 8 7]));

            qt = quatmultiply(obj.NomState(7:10)',dq')';

            dx_noRot   = [dx(1:6)           ; dx(10:15)];
            Xnom_noRot = [obj.NomState(1:6) ; obj.NomState(11:16)]; 

            CorrectedState_noRot = Xnom_noRot + dx_noRot;  
            CorrectedState       = [CorrectedState_noRot(1:6) ; qt ; CorrectedState_noRot(7:12)];

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
            % Updates the IMU and GPS readings
            % GRAVITYVECTOR = [0; 0; -9.81];  % ENU
            % rng(1,"twister");


            obj.KINPerIMU = obj.KINPerIMU + 1;


            % Get the true position, velocity, and acceleration in
            % coordinate frame.
            trueTranslationalVelocity = trueVelocity(1:3);
            trueAngularVelocity = [0; 0; trueVelocity(4)]; % the order of angular velocity should be X,Y,Z
            % trueAngularVelocity = [trueVelocity(4) ; 0 ; 0]; % the order of angular velocity should be X,Y,Z

            trueAcceleration = (trueTranslationalVelocity - obj.prevTrueVelocity(1:3))/obj.kindt;
            obj.prevTrueVelocity = trueVelocity;

            % Update IMU reading: coordinateAcceleration = GRAVITYVECTOR' - accelReading
            obj.KINRateIMU = obj.kinFs/obj.imuFs;
            if mod(obj.KINPerIMU,obj.KINRateIMU) == 0        
                [accelReading_body, gyroReading_body] = obj.IMU(trueAcceleration',trueAngularVelocity',obj.nedROTbody);
                obj.IMUreading  = {accelReading_body, gyroReading_body};

                [accelReading_body, gyroReading_body] = obj.IMUTruth(trueAcceleration',trueAngularVelocity',obj.nedROTbody);
                obj.IMUTruthreading  = {accelReading_body, gyroReading_body};
    
                obj.moveINS

            end
        end



    end

end