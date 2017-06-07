% Test the functionality of the PrincipalAxesOpt class

addpath('../');

% Select test case
% testCase = 'TA';
testCase = 'AS1';

% Test parameters - use similar inertias to SPHERES, but starting at t0=0
if strcmp(testCase,'TA')
    J = [1.239 1.1905]';                                        % tri-axial
    omegaB0_B = [0.939392242898362   0   0.500486277097766]';   % tri-axial
elseif strcmp(testCase,'AS1')
    J = [1.8534 1]';                                        	% axis-symmetric 1
    omegaB0_B = [1 0.4 0.4]';                                   % axis-symmetric 1
end
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
% omegaB_B = rigidBodyRotation.predictOmega(t);                          % noiseless
% omegaB_B = rigidBodyRotation.predictOmega(t) + normrnd(0,0.04,3,nT);    % noisy

% Create a rotation matrix to misalign the data
thetaTest = -[1 1 1]'*pi/4;
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

% Plot true data, and axis alignment data
figure(1); clf;
subplot(3,2,[1 3 5]);
supressB = 1;    % whether or not to supress assigned body axes
principalAxesOpt.plot(supressB,0);
plot3(omegaB_Bt(1,:),omegaB_Bt(2,:),omegaB_Bt(3,:),'-c','LineWidth',2);
RBtoBhat = RBtoGe' * RBtoG;
plot3([0 RBtoBhat(1,1)],[0 RBtoBhat(2,1)],[0 RBtoBhat(3,1)],'r');
plot3([0 RBtoBhat(1,2)],[0 RBtoBhat(2,2)],[0 RBtoBhat(3,2)],'g');
plot3([0 RBtoBhat(1,3)],[0 RBtoBhat(2,3)],[0 RBtoBhat(3,3)],'b');
% view(-75,45);
if ~supressB
    legend('Measured \omega_B','Estimated x_E','Estimated y_E','Estimated z_E', ...
        'Estimated x_B','Estimated y_B','Estimated z_B','Initial Conic Fit 1',...
        'Initial Conic Fit 2','Initial Conic Fit 3','Final Conic Fit 1','Final Conic Fit 2',...
        'Final Conic Fit 3','True \omega_B','True x_B','True y_B','True z_B', ...
        'Location','eastoutside'); 
    title(['Principal Axes Alignment. Error = ' num2str(rad2deg(rotErr)) '^o.']);
else
    legend('Measured \omega_B','Estimated x_E','Estimated y_E','Estimated z_E', ...
        'Inital Conic Fit 1','Initial Conic Fit 2','Initial Conic Fit 3',...
        'True \omega_B','True x_B','True y_B','True z_B','Location','eastoutside');
    title('Principal Axes Alignment');
end
xlabel('{}^B\omega_1 [rad/s]'); ylabel('{}^B\omega_2 [rad/s]'); zlabel('{}^B\omega_3 [rad/s]');

% Plot angular velocities projected by the proposed alignment
subplot(3,2,2);
plot(t,omegaB_Bt(1,:),'-c','LineWidth',2); hold on; grid on;
plot(t,omegaB_Be(1,:),'Color',[0.5 0 0.5],'LineWidth',1);
legend('True','Aligned Meas.');
title('Measurements Aligned Using Estimated {}_B^GR'); ylabel('{}^B\omega_1 [rad/s]');
subplot(3,2,4);
plot(t,omegaB_Bt(2,:),'-c','LineWidth',2); hold on; grid on;
plot(t,omegaB_Be(2,:),'Color',[0.5 0 0.5],'LineWidth',1); ylabel('{}^B\omega_2 [rad/s]');
subplot(3,2,6);
plot(t,omegaB_Bt(3,:),'-c','LineWidth',2); hold on; grid on;
plot(t,omegaB_Be(3,:),'Color',[0.5 0 0.5],'LineWidth',1); ylabel('{}^B\omega_3 [rad/s]');
xlabel('Time t [s]');