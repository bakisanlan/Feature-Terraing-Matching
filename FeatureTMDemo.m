clc;clear;close all;
rng(32,"twister");

IMUtype = 2;
IMUTs = 0.1; 
KINTs = 0.1;
initialState = zeros(16,1); refLLA = [41 29 0];
ReferenceFrame = 'ENU';
hAircraft = INSbot(initialState, IMUTs, KINTs, ReferenceFrame, IMUtype);

%% DTED
area = 'ITU';
hAIM                    = AerialImageModel(area);
hUAVScanner             = DatabaseScanner;
hReferenceMapScanner    = DatabaseScanner;

%% Scanner Settings

hUAVScanner.AIM          = hAIM;
hReferenceMapScanner.AIM = hAIM;

%% Aircraft Settings
x0      = 100;
y0      = 100;
z0      = 100;
pos0    = [x0 ; y0 ; z0];
psi0    = -20*pi/180;  %rad
aircraft_model_flag = 1; % Kinematic model
u = repmat([10 -2*pi/500],1000000,1);
% u = repmat([100 0*pi/500],10000,1);

% V0 = [100 0 0];
V0 = [u(1,1)*cos(psi0) ; u(1,1)*sin(psi0) ; 0];

hAircraft.TrueKinState          = [x0; y0; z0; psi0; 0; 0];
hAircraft.TrueKinQuat          = eul2quat(hAircraft.TrueKinState(4:6)')';
quat0 = eul2quat(hAircraft.TrueKinState(4:6)');

%% Set INSFILTER PROPORTIES
AccelerometerBias =	hAircraft.IMU.Accelerometer.ConstantBias;
GyroscopeBias     = hAircraft.IMU.Gyroscope.ConstantBias;

hAircraft.NomState  = [pos0 ; V0 ; quat0' ; AccelerometerBias' ; GyroscopeBias'];
hAircraft.TrueState = [pos0 ; V0 ; quat0' ; zeros(3,1)         ; zeros(3,1)];

hAircraft.prevTrueVelocity = [V0 ; u(1,2)];

%% Estimator settings
N_part = 1;                % number of particles
mu_part    = zeros(3,1);
std_part   = [10 ; 10 ; 0.1];
% std_part   = [0 ; 0 ; 0];

mu_kalman  = zeros(12,1);
cov_kalman = diag([std_part.^2 ; ones(9,1)]);
% exp_rate = 0;               % exploration rate of Particle Filter
R = 9;                % radar altitude standart error in meters
rayCast = false;
hEstimator = StateEstimatorMPF(N_part,mu_part,std_part,mu_kalman,cov_kalman,IMUTs,R);
rng(32,"twister");
hEstimator.hReferenceMapScanner = hReferenceMapScanner;
hEstimator.Accelerometer = hAircraft.IMU.Accelerometer;
hEstimator.Gyroscope     = hAircraft.IMU.Gyroscope;

%% Create trace arrays
tracePose            = [];
traceState           = [];
traceINSnoFusePose   = [];
traceINSPredict      = [];
tracePFPose          = [];
traceEstimatedPose   = [];

traceRotEul          = [];
traceEstimatedRotEul = [];
traceEstimatedState  = [];

particles_history(1:N_part,:) = hEstimator.particles(1:3,:)';

simTime = 0;
Tf = 80;
%Kinematic Model Simulation
inputParticle = cell(1);
AircraftNomState = cell(1);
UAVImage = cell(1);
true_rot = cell(1);
true_pos = cell(1);

j = 1;
ndownsample = 1; % for particle location to radar point location plot

while simTime < Tf
    idx = round(simTime/IMUTs + 1);
    hAircraft.move(u(idx,:),aircraft_model_flag);

    if mod(hAircraft.KINPerIMU,hAircraft.KINRateIMU) == 0

        AircraftNomState{j} = hAircraft.NomState;

        inputParticle{j} =[hAircraft.IMUreading{1}' ; hAircraft.IMUreading{2}'];

        % UAV camera scan
        UAVImage{j} = hUAVScanner.snapUAVImage(hAircraft.TrueKinState(1:2),rad2deg(hAircraft.TrueKinState(4)));

        true_rot{j} = [hAircraft.TrueKinState(4) 0 0];
        true_pos{j} = [hAircraft.TrueKinState(1:3)];

        truePose      = hAircraft.TrueKinState(1:3)';
        trueState     = hAircraft.TrueKinState';
        INSnoFusePose = hAircraft.NomState(1:3)';

        tracePose            = [tracePose            ; truePose];
        traceState           = [traceState           ; trueState];
        traceINSnoFusePose   = [traceINSnoFusePose   ; INSnoFusePose];

        j = j + 1;        
    end

    simTime = simTime + KINTs;

end
count_i = 0;

%Estimation Loop
for i = 1:j-1
    

    tic;
    param = hEstimator.getEstimate(inputParticle{i},AircraftNomState{i},UAVImage{i});
    toc;

    hAircraft.NomState = AircraftNomState{i};
    estState = hAircraft.correctINS(param.State);
    estState = [estState(1:6) ; rad2deg(quat2eul(estState(7:10)'))' ; estState(11:16)];
    estPose  = estState(1:3)';

    traceEstimatedPose   = [traceEstimatedPose   ; estPose];
    traceEstimatedState  = [traceEstimatedState  ; estState'];

    if mod((i-1),ndownsample) == 0

        particles_history(1+N_part*(count_i):N_part*(count_i+1),:) = AircraftNomState{i}(1:3)'+ hEstimator.particles(1:3,:)';
        count_i = count_i +1;
    end
    %%
    
    % pause(0.01)% particle_hist = hAircraft.NomState(1:3)'+ hEstimator.particles(1:3,:)';
    % figure;
    % clf;
    % p1 = plot3(particle_hist(:,1),particle_hist(:,2),particle_hist(:,3),'.k',MarkerSize=0.5);
    % p1.Color = [0.48,0.47,0.47];
    % hold on;
    % plot3(truePose(:,1),truePose(:,2),truePose(:,3),'*g')
    % plot3(INSnoFusePose(:,1),INSnoFusePose(:,2),INSnoFusePose(:,3),'^c')
    % plot3(estPose(:,1),estPose(:,2),estPose(:,3),'xr')
    % p.Color = [0.9290 0.6940 0.1250];
    % 
    % daspect([1 1 0.05]);
    % pbaspect([1 1 1]);
    % xlabel('x')
    % ylabel('y')
    % zlabel('h','Rotation',0)
    % lim = 200;
    % xlim([truePose(1) - lim truePose(1) + lim])
    % ylim([truePose(2) - lim truePose(2) + lim])
    % % legend('Particles','Ground Truth','INS without GPS/TerrainMatch','INS Predict(before TerrainMatch update)','PF Estimate','INS Measurement Update','Location','northeast')
    % grid on
    % view(2)
   
end


%% Error Calculation
valid_index = ~isnan(tracePose);

errorsINS   = vecnorm(tracePose - traceINSnoFusePose,2,2);
errorsFused = vecnorm(tracePose - traceEstimatedPose,2,2);
errorsFusedRot = vecnorm(rad2deg(traceState(:,4:6)) - (traceEstimatedState(:,7:9)),2,2);

mean_errorINS   = mean(errorsINS);
mean_errorFused = mean(errorsFused);
mean_errorFusedRot = mean(errorsFusedRot);

mean_stdINS          = sqrt(sum((errorsINS - mean_errorINS).^2)/N_part);
mean_stdFused        = sqrt(sum((errorsFused - mean_errorFused).^2)/N_part);
mean_stdFusedRot        = sqrt(sum((errorsFusedRot - mean_errorFusedRot).^2)/N_part);

disp(['Estimation error of INS ',num2str(mean_errorINS),' meters mean and ',num2str(mean_stdINS),' std'])
disp(['Estimation error of Fused Filter ',num2str(mean_errorFused),' meters mean and ',num2str(mean_stdFused),' std'])
disp(['Estimation error Rotation of Fused Filter ',num2str(mean_errorFusedRot),' degrees mean and ',num2str(mean_stdFusedRot),' std'])

%%
particles_history(:,2) = -particles_history(:,2) + 200;
particles_history_px   =  round(particles_history * (1/hAIM.mp));

tracePose(:,2) = -tracePose(:,2) + 200;
tracePose_px   =  round(tracePose * (1/hAIM.mp));

traceINSnoFusePose(:,2) = -traceINSnoFusePose(:,2) + 200;
traceINSnoFusePose_px   =  round(traceINSnoFusePose * (1/hAIM.mp));

traceEstimatedPose(:,2) = -traceEstimatedPose(:,2) + 200;
traceEstimatedPose_px   =  round(traceEstimatedPose * (1/hAIM.mp));

%%
figure
imshow(hAIM.I)
hold on;
% 
p1 = plot(particles_history_px(:,1),particles_history_px(:,2),'.y',MarkerSize=5);
% p1.Color = [0.48,0.47,0.47];
plot(tracePose_px(1:ndownsample:end,1),tracePose_px(1:ndownsample:end,2),'*g')
plot(traceINSnoFusePose_px(1:ndownsample:end,1),traceINSnoFusePose_px(1:ndownsample:end,2),'^c')
plot(traceEstimatedPose_px(1:ndownsample:end,1),traceEstimatedPose_px(1:ndownsample:end,2),'xr')


% p1 = plot(particles_history(:,1),particles_history(:,2),'.k',MarkerSize=0.5);
% p1.Color = [0.48,0.47,0.47];
% plot(tracePose(1:ndownsample:end,1),tracePose(1:ndownsample:end,2),'*g')
% plot(traceINSnoFusePose(1:ndownsample:end,1),traceINSnoFusePose(1:ndownsample:end,2),'^c')
% plot(traceEstimatedPose(1:ndownsample:end,1),traceEstimatedPose(1:ndownsample:end,2),'xr')


daspect([1 1 0.05]);
pbaspect([1 1 1]);
xlabel('x')
ylabel('y')
legend('Particles','Ground Truth','INS without GPS/TerrainMatch','PF Estimate','Sea Region','Normalized Terrain Uniqueness','Location','northeast')
grid on
view(2)
