% Test the functionality of the InertiaConstraints class

addpath('../');

% Select test case
% testCase = 'TA';
testCase = 'AS1';
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

%{

% Test parameters - use similar inertias to SPHERES, but starting at t0=0
J = [1.239 1.1905]';
J1 = J(1); J2 = J(2); J3 = 1;
RB0toW = eye(3);
omegaB0_B = [0.939392242898362   0   0.500486277097766]';
omegaB0_W = RB0toW * omegaB0_B;

% Create a test canonical rigid body
RigidBodyRotation(J,RB0toW,omegaB0_B,'omega0',0);
Ek = rigidBodyRotation.Ek; 
h = norm(rigidBodyRotation.getAngularMomentum());
%}

% Quarter-period and times [s]
T = rigidBodyRotation.T;
dt = 0.5;
t = 0:dt:2*T;
nT = length(t);

% Get the angular velocities, and rotation matrices
% omegaB_B = rigidBodyRotation.predictOmega(t);                          % noiseless
omegaB_B = rigidBodyRotation.predictOmega(t) + normrnd(0,0.04,3,nT);    % noisy

% Create the inertia constraints object
inertiaConstraints = InertiaConstraints(omegaB_B,rigidBodyRotation.energyState,....
                        rigidBodyRotation.inertiaSymmetry);

% Plot the inertia ratio constraints
figure(2); clf;
inertiaConstraints.plot(); grid on; axis equal;
plot(J(1),J(2),'*k','MarkerSize',10,'LineWidth',2);
xlim([0 4]); ylim([0 4]); xlabel('J_1'); ylabel('J_2');
title('Constraints on Inertia Ratios');

% Find the closest point for each equality constraint
Jtest = [1.5 2]';
% plot(Jtest(1),Jtest(2),'or');
for i = 1:6
    [pOpt,Jopt] = inertiaConstraints.findClosestPoint(Jtest,i);
%     plot(Jopt(1),Jopt(2),'oc');
%     plot([Jtest(1) Jopt(1)],[Jtest(2) Jopt(2)],'-r');
end

if strcmp(testCase,'TA')
    legend('Feasible Inertia Ratios','Convention Inequality Constraint J_1 \geq J_2',...
        'Convention Inequality Constraint J_2 \geq J_3 = 1', ...
        'Triangle Inequality Constraint J_2 + J_3 \geq J_1','High Energy / Low Energy Divider', ...
        '\alpha_1 Conic Equality Constraint', ... %'\gamma_1 Conic Equality Constraint', ...
        '\alpha_2 Conic Equality Constraint', ... %'\gamma_2 Conic Equality Constraint', ...
        '\alpha_3 Conic Equality Constraint', ... %'\gamma_3 Conic Equality Constraint 6', ...
        'Optimal Inertia Ratios','Location','NorthWest');
elseif strcmp(testCase,'AS1') || strcmp(testCase,'AS3')
    legend('Feasible Inertia Ratios','Convention Inequality Constraint J_1 \geq J_2',...
        'Convention Inequality Constraint J_2 \geq J_3 = 1', ...
        'Triangle Inequality Constraint J_2 + J_3 \geq J_1','High Energy / Low Energy Divider', ...
        'Linear Equality Constraint', 'Optimal Inertia Ratios','Location','NorthWest');
end
%{

legend('Linear Inequality Constraint 1','Linear Inequality Constraint 2', ...
    'Linear Inequality Constraint 3','Energy Conic Inequality Constraint', ...
    'Ang Vel Conic Equality Constraint 1','Ang Vel Conic Equality Constraint 2', ...
    'Ang Vel Conic Equality Constraint 3','Ang Vel Conic Equality Constraint 4', ...
    'Ang Vel Conic Equality Constraint 5','Ang Vel Conic Equality Constraint 6', ...
    'Actual Inertia Ratios','Test Point');
%}