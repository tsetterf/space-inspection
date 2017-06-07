% This script aims to do the basic and extended tests of the RigidBodyRotation class to ensure that 
% it is working. 

% It does this by getting the angular velocity solutions analytically, integrating them 
% numerically, and comparing the integrated orientation results with the analytic orientation
% results.

% Add path to rigid body dynamics
addpath('../');

%% Tri-Axial Test Case ===============

% Test parameters
ratiosSph = [1.239 1.1905]';               % approx fit SPHERES inertia ratios
% ratiosSph = [1.1990 1.1642]';               % SPHERES inertia ratios for OL_INT test
RB0toW = eye(3);                            % initial orientation

% Create LE case by choosing similar velocities to maximal during a SPHERES experiment
disp('Creating rigid body rotations ==============');
disp('LE RigidBodyRotation ---');
% omegaB0_BLE = [0.052359877559830 1.047197551196598 0.052359877559830]'; % from OL_INT test TS63
omegaB0_BLE = [0.939392242898362   0   0.500486277097766]';
rbrLE        = RigidBodyRotation(ratiosSph,RB0toW,omegaB0_BLE,'omega0',0);
% Create ME case by using the fact that Ek = h^2/(2*J2)
disp('ME RigidBodyRotation ---');
rbrME        = RigidBodyRotation(ratiosSph,RB0toW, ...
                [rbrLE.h^2/(2*ratiosSph(2)) rbrLE.h 0]','Ekht0',0);
% Create HE case by rearranging LE angular velocities 
disp('HE RigidBodyRotation ---');
omegaB0_BHE = [0 0.939392242898362 0.500486277097766]';
rbrHE        = RigidBodyRotation(ratiosSph,RB0toW,omegaB0_BHE,'omega0',0);
% Create LE case with principal rotation by specifying an angular velocity about the maj axis
disp('LE Principal RigidBodyRotation ---');
omegaB0_BLEP = [norm(omegaB0_BLE) 0 0]';
rbrLEP       = RigidBodyRotation(ratiosSph,RB0toW,omegaB0_BLEP,'omega0',0);

% Choose which case to use ***
rigidBodyRotation = rbrLE;
% rigidBodyRotation = rbrME;
% rigidBodyRotation = rbrHE;
% rigidBodyRotation = rbrLEP;

% Times
nT = 1000;
t = linspace(0,4*rigidBodyRotation.T,nT);

% Setup variable histories         
qWtoB = zeros(4,nT);     
qWtoB(:,1) = rot2quat(RB0toW);
orientErr = zeros(1,nT);

% Get the angular velocities, rotation matrices, and quaternions
omegaB_B = rigidBodyRotation.predictOmega(t);
RBtoW    = rigidBodyRotation.predictOrientation(t);
qWtoBint = integrateOmega(t,omegaB_B);
qErr = zeros(4,nT); qErr(:,1) = [0 0 0 1]';
for i = 2:nT
    qWtoB(:,i) = rot2quat(RBtoW(:,:,i)');
    qErr(:,i) = quatmult(qWtoB(:,i),quatconj(qWtoBint(:,i)));
    orientErr(i) = abs(rad2deg(wrapToPi(2*acos(min(1,max(-1,qErr(4,i)))))));
end

% Fix quaternions
qWtoB = quatfix(qWtoB);
qErr = quatfix(qErr);

% Perform self-consistency check to make sure that createMode 'Ekht0' works as well
disp('Performing self-consistency check =========');
Ek = rigidBodyRotation.Ek; h = rigidBodyRotation.h; t0 = rigidBodyRotation.t0;
rigidBodyRotation2 = RigidBodyRotation(ratiosSph,RB0toW,[Ek h t0]','Ekht0',0);
omegaSelfCons = sum(sum(omegaB_B - rigidBodyRotation2.predictOmega(t)));
orientationSelfCons = sum(sum(sum(RBtoW - rigidBodyRotation2.predictOrientation(t))));
if omegaSelfCons < 10e-6 && orientationSelfCons < 10e-6
    disp('Self-consistency check passed!')
else
    disp('Self-consistency check failed!')
    disp(['omegaSelfCons = ' num2str(omegaSelfCons)]);
    disp(['orientationSelfCons = ' num2str(orientationSelfCons)]);
end

% Plot the angular velocities
figure(1); clf;
subplot(3,1,1);
plot(t,omegaB_B(1,:),'r'); hold on; grid on;
ylabel('\omega_1 [rad/s]'); axis([0 max(t) -2 2]);
title('RigidBody Analytic Angular Velocities');
subplot(3,1,2);
plot(t,omegaB_B(2,:),'r'); hold on; grid on;
ylabel('\omega_2 [rad/s]'); axis([0 max(t) -2 2]);
subplot(3,1,3);
plot(t,omegaB_B(3,:),'r'); hold on; grid on;
ylabel('\omega_3 [rad/s]'); xlabel('t'); axis([0 max(t) -2 2]);

% Plot the quaternions
figure(2); clf;
subplot(5,1,1);
plot(t,qWtoB(1,:),'r'); hold on; grid on;
plot(t,qWtoBint(1,:),'--b'); ylabel('q_1'); axis([0 max(t) -1 1]); 
title('RigidBody Analytic Orientation and Integrated Orientation');
legend('Analytic','Integrated from analytic {}^B\omega');
subplot(5,1,2);
plot(t,qWtoB(2,:),'r'); hold on; grid on;
plot(t,qWtoBint(2,:),'--b'); ylabel('q_2'); axis([0 max(t) -1 1]);
subplot(5,1,3);
plot(t,qWtoB(3,:),'r'); hold on; grid on;
plot(t,qWtoBint(3,:),'--b'); ylabel('q_3'); axis([0 max(t) -1 1]);
subplot(5,1,4);
plot(t,qWtoB(4,:),'r'); hold on; grid on;
plot(t,qWtoBint(4,:),'--b'); ylabel('q_4'); axis([0 max(t) -1 1]);
subplot(5,1,5);
plot(t,orientErr,'-m'); hold on; grid on; 
ylabel('\theta_{diff} [deg]');

disp(['Average quaternion error ' num2str(mean(orientErr)) ' deg']);

% Calculate position of angular momentum in inertial frame
h_W = zeros(3,nT); hint_W = zeros(3,nT);
for i = 1:nT
    h_W(:,i) = quat2rot(qWtoB(:,i))' * diag(rigidBodyRotation.J) * omegaB_B(:,i);
    hint_W(:,i) = quat2rot(qWtoBint(:,i))' * diag(rigidBodyRotation.J) * omegaB_B(:,i);
end
hact_W = rigidBodyRotation.getAngularMomentumVector();

% Plot a check of the angular momentum
figure(3); clf;
subplot(3,1,1);
plot(t,h_W(1,:),'r'); hold on; grid on;
plot(t,hint_W(1,:),'--b'); plot([t(1) t(end)],[hact_W(1) hact_W(1)],'--g'); 
ylabel('{}^Wh_1'); title('Check for constant angular momentum vector'); 
legend('Analytic','Integrated from {}^B\omega','Actual Assigned');
subplot(3,1,2);
plot(t,h_W(2,:),'r'); hold on; grid on;
plot(t,hint_W(2,:),'--b'); plot([t(1) t(end)],[hact_W(2) hact_W(2)],'--g'); 
ylabel('{}^Wh_2');
subplot(3,1,3);
plot(t,h_W(3,:),'r'); hold on; grid on;
plot(t,hint_W(3,:),'--b'); plot([t(1) t(end)],[hact_W(3) hact_W(3)],'--g'); 
ylabel('{}^Wh_3');

% Plot the polhode and ellipsoids
figure(4); clf;
rigidBodyRotation.plotPolhode(1,1);
legend('Momentum Ellipsoid','Energy Ellipsoid','Polhode','Inital Angular Velocity','Location','NorthWest'); 
title('Tri-Axial Polhode');

