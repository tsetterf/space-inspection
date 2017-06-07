% Test the complete process of polhode estimation using PrincipalAxesOpt and InertiaRatiosOpt

addpath('../');

% Select test case
testCase = 'TA'; 
% testCase = 'AS1';
% testCase = 'AS3';

% Test parameters - use similar inertias to SPHERES, but starting at t0=0
if strcmp(testCase,'TA')
    J = [1.239 1.1905]';                                        % tri-axial
    omegaB0_B = [0.939392242898362   0   0.500486277097766]';   % tri-axial
elseif strcmp(testCase,'AS1')
    J = [1.8534 1]';                                        	% axis-symmetric 1
    omegaB0_B = [1 0.4 0.4]';                                   % axis-symmetric 1
elseif strcmp(testCase,'AS3')
    J = [1.6839 1.6839]';                                       % axis-symmetric 3
    omegaB0_B = [0.4 0.4 1]';                                   % axis-symmetric 3
end

%% Principal Axes Alignment

% Test parameters - use similar inertias to SPHERES, but starting at t0=0
J1 = J(1); J2 = J(2); J3 = 1;
RB0toW = eye(3);

% Create a test canonical rigid body
rigidBodyRotation = RigidBodyRotation(J,RB0toW,omegaB0_B,'omega0',0);
Ek = rigidBodyRotation.Ek; 
h = norm(rigidBodyRotation.getAngularMomentum());

% Quarter-period and times [s]
T = rigidBodyRotation.T;
dt = 0.5;
t = 0:dt:2.5*T;
tExt = 0:dt:4*t(end);
nT = length(t);

% Get the angular velocities for the sampling time, and extended (for truth)
omegaB_Bt = rigidBodyRotation.predictOmega(t);                          % noiseless
omegaB_BtExt = rigidBodyRotation.predictOmega(tExt);                    % noiseless, extended
% omegaB_B = rigidBodyRotation.predictOmega(t);                          % noiseless
omegaB_B = rigidBodyRotation.predictOmega(t) + normrnd(0,0.04,3,nT);    % noisy
RBtoW    = rigidBodyRotation.predictOrientation(t);

% Create a rotation matrix to misalign the data
thetaTest = [1 1 1]'*pi/4;
RBtoG = Exp(thetaTest);
RGtoB = RBtoG';

% Get the rotated angular velocities and rotation from geometric to world frame
omegaB_Gt = RGtoB' * omegaB_Bt;
omegaB_G  = RGtoB' * omegaB_B;

% Create the principal axes optimizer
principalAxesOpt = PrincipalAxesOpt(omegaB_G);

% Perform optimization
[RBtoGe,costOpte] = principalAxesOpt.optimize();

% Get magnitude of rotation error
rotErr = norm(Log( RBtoG'*RBtoGe ));

% Get the estimated aligned angular velocities
omegaB_Be = principalAxesOpt.RBtoG' * omegaB_G;

% Get the estimated orientations
RBtoWe = zeros(size(RBtoW));
for i = 1:size(RBtoW,3)
    RGitoW = RBtoW(:,:,i) * RGtoB;
    RBtoWe(:,:,i) = RGitoW * RBtoGe;
end

% Plot true data, and axis alignment data
figure(1); clf;
subplot(3,2,[1 3 5]);
principalAxesOpt.plot(0,1);
plot3(omegaB_Gt(1,:),omegaB_Gt(2,:),omegaB_Gt(3,:),'-k','LineWidth',2);
plot3([0 RBtoG(1,1)],[0 RBtoG(2,1)],[0 RBtoG(3,1)],'r');
plot3([0 RBtoG(1,2)],[0 RBtoG(2,2)],[0 RBtoG(3,2)],'g');
plot3([0 RBtoG(1,3)],[0 RBtoG(2,3)],[0 RBtoG(3,3)],'b');
view(-75,45);
legend('Measured {}^G\omega_B','Estimated {}^Gx_E','Estimated {}^Gy_E','Estimated {}^Gz_E', ...
    'Estimated {}^Gx_B','Estimated {}^Gy_B','Estimated {}^Gz_B','Fit Conic 1','Fit Conic 2', ...
    'Fit Conic 3','True {}^G\omega_B','True {}^Gx_B','True {}^Gy_B','True {}^Gz_B', ...
    'Location','eastoutside');
xlabel('{}^G\omega_1 [rad/s]'); ylabel('{}^G\omega_2 [rad/s]'); zlabel('{}^G\omega_3 [rad/s]'); 
title(['Alignment of Principal Axes. Error = ' num2str(rad2deg(rotErr)) '^o']);

% Plot angular velocities projected by the proposed alignment
subplot(3,2,2);
plot(tExt,omegaB_BtExt(1,:),'-k','LineWidth',2); hold on; grid on;
plot(t,omegaB_Be(1,:),'Color',[0.5 0 0.5],'LineWidth',1); ylabel('{}^B\omega_1 [rad/s]');
subplot(3,2,4);
plot(tExt,omegaB_BtExt(2,:),'-k','LineWidth',2); hold on; grid on;
plot(t,omegaB_Be(2,:),'Color',[0.5 0 0.5],'LineWidth',1); ylabel('{}^B\omega_2 [rad/s]');
subplot(3,2,6);
plot(tExt,omegaB_BtExt(3,:),'-k','LineWidth',2); hold on; grid on;
plot(t,omegaB_Be(3,:),'Color',[0.5 0 0.5],'LineWidth',1); ylabel('{}^B\omega_3 [rad/s]');
xlabel('Time t [s]');

%% Inertia Ratios Optimization

% Use identity covariances to solve least squares problem
covOmegaB = repmat(eye(3),1,1,nT);

% Create the inertia ratios optimizer
inertiaRatiosOpt = InertiaRatiosOpt(omegaB_Be,covOmegaB,RBtoWe,t,principalAxesOpt.energyState,...
                    principalAxesOpt.inertiaSymmetry);

% Test the optimization
[pOpt,costOpt,indOpt] = inertiaRatiosOpt.optimize();
disp(['Actual J: ' num2str(J')]);

% Get the intialization angular velocities and the optimal ones (both from the best constraint)
% omegaB_BinitExt = inertiaRatiosOpt.predictOmega( ...
%                     inertiaRatiosOpt.p0(indOpt),inertiaRatiosOpt.t00(indOpt),tExt,indOpt);
% omegaB_BoptExt = inertiaRatiosOpt.predictOmega( ...
%                     inertiaRatiosOpt.pOpt(indOpt),inertiaRatiosOpt.t0opt(indOpt),tExt,indOpt);
omegaB_BinitExt = inertiaRatiosOpt.predictOmegaInit(tExt);
omegaB_BoptExt = inertiaRatiosOpt.predictOmegaOpt(tExt);
                
% Plot the optimized angular velocity profiles
figure(1);
subplot(3,2,2);
plot(tExt,omegaB_BinitExt(1,:),'--r'); plot(tExt,omegaB_BoptExt(1,:),'--g','LineWidth',2);
legend('True','Aligned Meas.','Init. Ang. Mom. Fit','Optimal Fit');
title('Fitting and Forward Prediction of Body Angular Velocities');
subplot(3,2,4);
plot(tExt,omegaB_BinitExt(2,:),'--r'); plot(tExt,omegaB_BoptExt(2,:),'--g','LineWidth',2);
subplot(3,2,6);
plot(tExt,omegaB_BinitExt(3,:),'--r'); plot(tExt,omegaB_BoptExt(3,:),'--g','LineWidth',2);