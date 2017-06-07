% Test the functionality of the InertiaRatioOpt class

addpath('../');

% Select test case
testCase = 'TA';
% testCase = 'AS1';
% testCase = 'AS3';

% Test parameters - use similar inertias to SPHERES, but starting at t0=0
if strcmp(testCase,'TA')
    J = [1.239 1.1905]';                                        % tri-axial
%     omegaB0_B = [0.3686 -0.9873 0.1224]';                       % different t0
    omegaB0_B = [0.939392242898362   0   0.500486277097766]';   % tri-axial
%     omegaB0_B = [0.3350 -1.0028 -0.0867]';                    % different t0
elseif strcmp(testCase,'AS1')
    J = [1.8534 1]';                                        	% axis-symmetric 1
    omegaB0_B = [1 0.4 0.4]';                                   % axis-symmetric 1
elseif strcmp(testCase,'AS3')
    J = [1.6839 1.6839]';                                       % axis-symmetric 3
    omegaB0_B = [0.4 0.4 1]';                                   % axis-symmetric 3
end

% Test parameters - use similar inertias to SPHERES, but starting at t0=0
J1 = J(1); J2 = J(2); J3 = 1;
RB0toW = eye(3);
omegaB0_W = RB0toW * omegaB0_B;

% Create a test canonical rigid body
rigidBodyRotation = RigidBodyRotation(J,RB0toW,omegaB0_B,'omega0',0);
Ek = rigidBodyRotation.Ek; 
h = norm(rigidBodyRotation.getAngularMomentum());

% Quarter-period and times [s]
T = rigidBodyRotation.T;
dt = 0.5;
t = 0:dt:2.5*T;
tExt = 0:dt:2*t(end);
nT = length(t);

% Get the angular velocities for the sampling time, and extended (for truth)
omegaB_Bt = rigidBodyRotation.predictOmega(t);                          % noiseless
omegaB_BtExt = rigidBodyRotation.predictOmega(tExt);                    % noiseless, extended
% omegaB_B = rigidBodyRotation.predictOmega(t);                         % noiseless
omegaB_B = rigidBodyRotation.predictOmega(t) + normrnd(0,0.04,3,nT);    % noisy
RBtoW    = rigidBodyRotation.predictOrientation(t);

% Use identity covariances to solve least squares problem
covOmegaB = repmat(eye(3),1,1,nT);

% Create the inertia ratios optimizer
inertiaRatiosOpt = InertiaRatiosOpt(omegaB_B,covOmegaB,RBtoW,t,rigidBodyRotation.energyState,...
                                        rigidBodyRotation.inertiaSymmetry);

% Test the optimization
disp('Testing optimization ...');
[pOpt,costOpt,indOpt] = inertiaRatiosOpt.optimize();
omegaB_BoptExt = zeros(3,length(tExt),3);

% Plot the cost vs parameter for the different inertia ratio constraints
figure(1); clf;
subplot(3,2,[1 3 5]);
inertiaRatiosOpt.plotCostVsP(1,'-m'); hold on; grid on;  
if strcmp(testCase,'TA')
    inertiaRatiosOpt.plotCostVsP(2,'--b');                 
    inertiaRatiosOpt.plotCostVsP(3,'-.c');     
end
xlabel('Parameter p [-]'); ylabel('Cost [rad/s]^2');

% Initialization omegas on optimal constraint
omegaB_Bopt0Ext = inertiaRatiosOpt.predictOmegaInit(tExt);

disp(['Actual J:  ' num2str(J')]);
disp(['Actual t0: ' num2str(rigidBodyRotation.t0)]);
    
% Optimized omegas
if strcmp(testCase,'TA')
    for i = 1:3
        omegaB_BoptExt(:,:,i) = inertiaRatiosOpt.predictOmega( ...
                                    inertiaRatiosOpt.pOpt(i),inertiaRatiosOpt.t0opt(i),tExt,i);
    end
    plot(inertiaRatiosOpt.p0(1),inertiaRatiosOpt.cost0(1),'om'); 
    plot(inertiaRatiosOpt.p0(2),inertiaRatiosOpt.cost0(2),'ob');
    plot(inertiaRatiosOpt.p0(3),inertiaRatiosOpt.cost0(3),'oc'); 
    plot(inertiaRatiosOpt.pOpt(1),inertiaRatiosOpt.costOpt(1),'*m'); 
    plot(inertiaRatiosOpt.pOpt(2),inertiaRatiosOpt.costOpt(2),'*b');
    plot(inertiaRatiosOpt.pOpt(3),inertiaRatiosOpt.costOpt(3),'*c'); 
    legend('Constraint \alpha_1','Constraint \alpha_2', ...
        'Constraint \alpha_3','p^0 Constraint \alpha_1','p^0 Constraint \alpha_2', ...
        'p^0 Constraint \alpha_3', 'p^* Constraint \alpha_1','p^* Constraint \alpha_2', ...
        'p^* Constraint \alpha_3');
title('Cost vs. Parameter p Along Equality Constraints.');

elseif strcmp(testCase,'AS1') || strcmp(testCase,'AS3')    
    omegaB_BoptExt(:,:,1) = inertiaRatiosOpt.predictOmega( ...
                                    inertiaRatiosOpt.pOpt(1),inertiaRatiosOpt.t0opt(1),tExt,1);
    plot(inertiaRatiosOpt.p0(1),inertiaRatiosOpt.cost0(1),'om'); 
    plot(inertiaRatiosOpt.pOpt(1),inertiaRatiosOpt.costOpt(1),'*m');
    legend('Constraint 1','p^0 Constraint 1','p^* Constraint 1');
	title('Cost vs. Parameter p Along Equality Constraint.');
end

% Plot the fits to angular velocity
figure(1);
subplot(3,2,2);
plot(tExt,omegaB_BtExt(1,:),'k','LineWidth',2); hold on; grid on;
plot(t,omegaB_B(1,:),'Color',[0.5 0 0.5]);
plot(tExt,omegaB_Bopt0Ext(1,:),'--r');
plot(tExt,omegaB_BoptExt(1,:,1),'-m'); 
ylabel('\omega_1 [rad/s]'); 
if strcmp(testCase,'TA')
    title(['Angular Velocity Fits. Min Cost: ' ...
            num2str(inertiaRatiosOpt.costOpt(inertiaRatiosOpt.indOpt)) ...
           ', Best Constraint: \alpha_' num2str(inertiaRatiosOpt.indOpt) '.']);
    plot(tExt,omegaB_BoptExt(1,:,2),'--b'); plot(tExt,omegaB_BoptExt(1,:,3),'-.c'); 
    legend('True','Measured','Initialization','Fit Constr. \alpha_1','Fit Constr. \alpha_2', ...
        'Fit Constr. \alpha_3');
elseif strcmp(testCase,'AS1') || strcmp(testCase,'AS3')
    title(['Angular Velocity Fits. Min Cost: ' ...
            num2str(inertiaRatiosOpt.costOpt(inertiaRatiosOpt.indOpt)) '.']);
    legend('True','Measured','Initialization','Fit Constr. 1');
end

subplot(3,2,4);
plot(tExt,omegaB_BtExt(2,:),'k','LineWidth',2); hold on; grid on;
plot(t,omegaB_B(2,:),'Color',[0.5 0 0.5]);
plot(tExt,omegaB_Bopt0Ext(2,:),'--r');
plot(tExt,omegaB_BoptExt(2,:,1),'-m'); 
ylabel('\omega_2 [rad/s]'); 
if strcmp(testCase,'TA')
    plot(tExt,omegaB_BoptExt(2,:,2),'--b'); plot(tExt,omegaB_BoptExt(2,:,3),'-.c'); 
end

subplot(3,2,6);
plot(tExt,omegaB_BtExt(3,:),'k','LineWidth',2); hold on; grid on;
plot(t,omegaB_B(3,:),'Color',[0.5 0 0.5]);
plot(tExt,omegaB_Bopt0Ext(3,:),'--r');
plot(tExt,omegaB_BoptExt(3,:,1),'-m'); 
xlabel('Time t [s]'); ylabel('\omega_3 [rad/s]');
if strcmp(testCase,'TA')
    plot(tExt,omegaB_BoptExt(3,:,2),'--b'); plot(tExt,omegaB_BoptExt(3,:,3),'-.c'); 
end

% Plot omegaPar vs parameter
% figure(2); clf;
% inertiaRatiosOpt.plotOmegaParVsP(1,'r'); grid on;
% ylabel('Parallel angular velocity magnitude ||\omega_{||}|| [rad/s]');
% xlabel('Parameter p [-]');
% title(['Parallel angular velocity magnitude \omega_{||} vs. Parameter p Along ' ...
%    'Various Equality Constraints']);
% if strcmp(testCase,'TA')
%     inertiaRatiosOpt.plotOmegaParVsP(3,'--b');
%     inertiaRatiosOpt.plotOmegaParVsP(5,'-.c'); 
% end

