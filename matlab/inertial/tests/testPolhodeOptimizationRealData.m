% Test the functionality of the InertiaRatioOpt and obtain ground truth inertia ratios data
% using gyro data from ISS tests

addpath('../../v0.6');
addpath('../');

% Select test data and set maneuver to start using data, a buffer of data to ignore, 
% the approx quarter period T, the energy state and inertia symmetry
dataDir = 'D:\PhD\vmshare\data\TS53T7R3\';
manStart = 4;
buffTime = 30;
Tapprox = 29.25;

% dataDir = 'D:\PhD\vmshare\data\TS38T2R1\'; % TweddlePhD 73.3s - 130.1s
% manStart = 3;
% buffTime = 31.39;
% Tapprox = 14.2;

%%{

% Load the test data
load([dataDir 'spheresData.mat']);

% Get the gyro data in the imu frame (geometric frame)
time = data{2}.BackTel.StdTime;                         % best for TS53T7R3
% time = data{2}.BackTel.StdTime - data{2}.TestsTime;     % best for TS38T2R1
% figure(14); clf; plot(time); hold on;
time = time(time >= data{2}.SOH.TestTime(find(data{2}.SOH.ManNum==manStart,1)) + buffTime*1e3);
plot(time);
gyro = data{2}.BackTel.StdState(11:13,data{2}.BackTel.StdTime >= time(1)); % TS53T7R3
% gyro = data{2}.BackTel.StdState(11:13,data{2}.BackTel.StdTime - data{2}.TestsTime >= time(1)); % TS38T2R1
gyroNoBias = zeros(size(gyro));
quat = data{2}.BackTel.StdState(7:10,data{2}.BackTel.StdTime >= time(1));  % TS53T7R3
% quat = data{2}.BackTel.StdState(7:10,data{2}.BackTel.StdTime - data{2}.TestsTime >= time(1)); % TS38T2R1
time = (time - time(1))/1000;

% Split into angular velocity segments of duration 4*T
numSplit = floor(time(end)/(4*Tapprox));
gyro = gyro(:,time <= 4*Tapprox*numSplit);
quat = quat(:,time <= 4*Tapprox*numSplit);
quat = quatfix(quat);
time = time(time <= 4*Tapprox*numSplit);

% Colors
col(1:3,1) = [0 0   1]';
col(1:3,2) = [1 0.6 0]';
col(1:3,3) = [0.5 0 0.5]';
col(1:3,4) = [0 1 0]';
col(1:3,5) = [0 1 1]';

% Plot the polhode and solve for each segment
figure(1); clf;
RBtoGopt = zeros(3,3,numSplit); costPrincipalAxesOpt = zeros(1,numSplit);
Jopt = zeros(2,numSplit); costInertiaRatiosOpt = zeros(1,numSplit);
principalAxesOpt = {};
inertiaRatiosOpt = {};
polhodeLegend = {};
for i = 1:numSplit
    
    % Get angular velocities and orientations
    iC = mod(i-1,5) + 1;
    t = time(time >= (i-1)*4*Tapprox & time <= i*4*Tapprox);    
    nT = length(t);
    omegaB_G = gyro(:,time >= (i-1)*4*Tapprox & time <= i*4*Tapprox);
    qWtoB = quat(:,time >= (i-1)*4*Tapprox & time <= i*4*Tapprox);
    RBtoW = zeros(3,3,nT);
    for j = 1:nT
        RBtoW(:,:,j) = quat2rot(qWtoB(:,j))';
    end
    
    % Find and remove gyroscope biases. First, refer all quaternions to their first rotation
    qB0toB = zeros(4,nT);
    for j = 2:nT
        qB0toB(:,j) = quatmult(qWtoB(:,j),quatconj(qWtoB(:,1)));
    end
    qB0toB(:,1) = [0 0 0 1]';
    qB0toB = quatfix(qB0toB);
    
    % Initial guess for bias
    bias0 = [0 0 0]';

    % Find the best estimate for bias to match the measured (unbiased) quaternions
    % Do this by matching the integrated angular velocities with the metrology quaternion data
    [bias,~] = lsqcurvefit(@(bg,tm)integrateOmega(tm,omegaB_G-repmat(bg,1,nT)),bias0,t,qB0toB);
    
    disp(['Found bias: ' num2str(bias')]);
    
    % Remove the bias from the angular velocities
    omegaB_G = omegaB_G - repmat(bias,1,nT);
    gyroNoBias(:,time >= (i-1)*4*Tapprox & time <= i*4*Tapprox) = omegaB_G;
    
    % Plot the polhode
    figure(1);
    subplot(3,2,[1 3]);
    plot3(omegaB_G(1,:),omegaB_G(2,:),omegaB_G(3,:),'-k','Color',col(:,iC)); grid on;  hold on; 
    axis equal; view(135,35);
    polhodeLegend{i} = ['Polhode ' num2str((i-1)*4*Tapprox) '-' num2str(i*4*Tapprox) ' s'];
    
    % Create the principal axes optimizer
    principalAxesOpt{i} = PrincipalAxesOpt(omegaB_G);

    % Perform optimization
    [RBtoGopt(:,:,i),costOpt] = principalAxesOpt{i}.optimize();
    costPrincipalAxesOpt(i) = costOpt;
    
    figure(2*i); clf;
    principalAxesOpt{i}.plot();
    plot3([0 1],[0 0],[0 0],'-.r'); plot3([0 0],[0 1],[0 0],'-.g'); plot3([0 0],[0 0],[0 1],'-.b');    
    drawnow;
    
    
    % Get the estimated aligned angular velocities
    omegaB_B = principalAxesOpt{i}.RBtoG' * omegaB_G;
    
    % Use identity covariances to solve least squares problem
    covOmegaB = repmat(eye(3),1,1,nT);

    % Create the inertia ratios optimizer
    inertiaRatiosOpt{i} = InertiaRatiosOpt(omegaB_B,covOmegaB,RBtoW,t-(i-1)*4*Tapprox,...
                            principalAxesOpt{i}.energyState,principalAxesOpt{i}.inertiaSymmetry);

    % Test the optimization
    disp('Testing optimization ...');
    [pOpt,costOpt,indOpt] = inertiaRatiosOpt{i}.optimize();
    Jopt(:,i) = inertiaRatiosOpt{i}.getInertiaRatios();
    costInertiaRatiosOpt(i) = costOpt;
    
end

% Get consensus value. see logbook #3 pp 150
sumTheta = 0; sumThetaDenom = 0; sumJ = zeros(2,1); sumJDenom = 0; thetaOpt = zeros(3,numSplit);
for i = 1:numSplit
    thetaOpt(:,i) = Log(principalAxesOpt{i}.RBtoG);
    sumTheta = sumTheta + Log(principalAxesOpt{i}.RBtoG)/principalAxesOpt{i}.costOpt;
    sumThetaDenom = sumThetaDenom + 1/principalAxesOpt{i}.costOpt;
    sumJ = sumJ + inertiaRatiosOpt{i}.getInertiaRatios()/inertiaRatiosOpt{i}.getCost();
    sumJDenom = sumJDenom + 1/inertiaRatiosOpt{i}.getCost();
end
thetaCons = 1/sumThetaDenom * sumTheta;
RBtoGcons = Exp(thetaCons);  
xB_G = RBtoGcons(:,1); yB_G = RBtoGcons(:,2); zB_G = RBtoGcons(:,3); 
Jcons = sumJ/sumJDenom;

% Plot consensus frame B
figure(1);
subplot(3,2,[1 3]);
plot3([0 xB_G(1)],[0 xB_G(2)],[0 xB_G(3)],'--r','LineWidth',3);
plot3([0 yB_G(1)],[0 yB_G(2)],[0 yB_G(3)],'--g','LineWidth',3);
plot3([0 zB_G(1)],[0 zB_G(2)],[0 zB_G(3)],'--b','LineWidth',3);
plot3([0 1],[0 0],[0 0],'-.r'); plot3([0 0],[0 1],[0 0],'-.g'); plot3([0 0],[0 0],[0 1],'-.b');
polhodeLegend{end+1} = '^Gx_B'; polhodeLegend{end+1} = '^Gy_B'; polhodeLegend{end+1} = '^Gz_B'; 
polhodeLegend{end+1} = '^Gx_G'; polhodeLegend{end+1} = '^Gy_G'; polhodeLegend{end+1} = '^Gz_G'; 
title('SPHERES Polhode'); legend(polhodeLegend,'Location','EastOutside');

disp(['Optimal RGtoB angle-axis vectors thetaOpt with std dev = ' num2str(std(thetaOpt(1,:))) ...
        ',' num2str(std(thetaOpt(2,:))) ',' num2str(std(thetaOpt(3,:)))]);
thetaOpt
disp('Cost for Optimal Rotation Matrices RGtoB');
costPrincipalAxesOpt
disp(['Consensus value of thetaOpt = ' num2str(thetaCons')]);
disp('Consensus RBtoG = ');
RBtoGcons
disp(['Optimal Inertia Ratios [J1 J2] with std dev = ' num2str(std(Jopt(1,:))) ...
        ',' num2str(std(Jopt(2,:)))]);
Jopt
disp('Cost for Optimal Inertia Ratios [J1 J2]');
costInertiaRatiosOpt
disp(['Consensus value for Jopt = ' num2str(Jcons')]);
%}

% Output table
axOpt = zeros(size(thetaOpt)); angOpt = zeros(numSplit,1);
axCons = thetaCons/norm(thetaCons); angCons = rad2deg(norm(thetaCons));
Iesl = [2.41e-2 -1.3e-4 -1.42e-4; -1.3e-4 2.34e-2 5.74e-5; -1.42e-4 5.74e-5 2.01e-2];
[vecEsl,valEsl] = eig(Iesl);
Iesl = flipud(diag(valEsl));
Jesl = Iesl/Iesl(3); Jesl = Jesl(1:2);
RBtoGesl = fliplr(vecEsl); RBtoGesl(:,2) = -RBtoGesl(:,2);
thetaOptEsl = Log(RBtoGesl); 
axOptEsl = thetaOptEsl/norm(thetaOptEsl); angOptEsl = rad2deg(norm(thetaOptEsl));
k1 = 0.0517; k2 = 0.0971; 
Jtwed = [exp(k1) 1 exp(-k2)]'; Jtwed = Jtwed/Jtwed(3); Jtwed = Jtwed(1:2);
RBtoGtwed = [0.987 0 0.158; 0 1 0; -0.158 0 0.987];
thetaOptTwed = Log(RBtoGtwed); 
axOptTwed = thetaOptTwed/norm(thetaOptTwed); angOptTwed = rad2deg(norm(thetaOptTwed));
disp('=======================================');
disp('& $\mathbf{a}^* \; [-]$ & $\theta^* \\; [^\\circ]$ & Cost$_{\theta}$ & $J_1^* \; [-]$ & $J_2^* \; [-]$ & Cost$_J$ \\');
disp('\hline');
for i = 1:numSplit
    axOpt(:,i) = thetaOpt(:,i)/norm(thetaOpt(:,i)); angOpt(i) = rad2deg(norm(thetaOpt(:,i)));
    fprintf(1,'Polhode %d-%d s & $\\left[ %.4f, \\; %.4f, \\; %.4f \\right]$ & %.4f & %.4f & %.4f & %.4f & %.4f \\\\\n', ...
        [(i-1)*4*Tapprox i*4*Tapprox axOpt(:,i)' angOpt(i) costPrincipalAxesOpt(i) Jopt(:,i)' ...
            costInertiaRatiosOpt(i)]);
end
fprintf(1,'Standard deviation & $\\left[ %.4f, \\; %.4f, \\; %.4f \\right]$ & %.4f & & %.4f & %.4f & \\\\\n', ...
    [std(axOpt,0,2)' std(angOpt) std(Jopt,0,2)']);
disp('\hline');
fprintf(1,'Consensus values & $\\left[ %.4f, \\; %.4f, \\; %.4f \\right]$ & %.4f & & %.4f & %.4f & \\\\\n', ...
    [axCons' angCons Jcons']);
disp('\hline');
fprintf(1,'Eslinger values~\\cite{Eslinger2013} & $\\left[ %.4f, \\; %.4f, \\; %.4f \\right]$ & %.4f & & %.4f & %.4f & \\\\\n', ...
    [axOptEsl' angOptEsl Jesl']);
disp('\hline');
fprintf(1,'Tweddle values~\\cite{TweddlePhD} & $\\left[ %.4f, \\; %.4f, \\; %.4f \\right]$ & %.4f & & %.4f & %.4f & \\\\\n', ...
    [axOptTwed' angOptTwed Jtwed']);

%%{
    
% Complete main plot
omegaB_B0m = RBtoGcons' * gyroNoBias(:,1);
Ek0 = 0.5*(Jcons(1)*omegaB_B0m(1)^2 + Jcons(2)*omegaB_B0m(2)^2 + 1*omegaB_B0m(3)^2);
h0 = sqrt((Jcons(1)*omegaB_B0m(1))^2 + (Jcons(2)*omegaB_B0m(2))^2 + (1*omegaB_B0m(3))^2);
for i = 1:numSplit
    
    iC = mod(i-1,5) + 1;
    t = time(time >= (i-1)*4*Tapprox & time <= i*4*Tapprox); 
    
    % Aligned measurements
    figure(1);
    omegaB_Bm = RBtoGcons' * gyroNoBias(:,time >= (i-1)*4*Tapprox & time <= i*4*Tapprox);
    subplot(3,2,2); plot(t,omegaB_Bm(1,:),'Color',col(:,iC)); grid on; hold on; 
    ylabel('\omega_1 [rad/s]');
    subplot(3,2,4); plot(t,omegaB_Bm(2,:),'Color',col(:,iC)); grid on; hold on; 
    ylabel('\omega_2 [rad/s]');
    subplot(3,2,6); plot(t,omegaB_Bm(3,:),'Color',col(:,iC)); grid on; hold on; 
    ylabel('\omega_3 [rad/s]'); 
    xlabel('Time [s]');
    
    % Kinetic energy and angular momentum
    Ek = 0.5*(Jcons(1)*omegaB_Bm(1,:).^2 + Jcons(2)*omegaB_Bm(2,:).^2 + 1*omegaB_Bm(3,:).^2);
    h = sqrt((Jcons(1)*omegaB_Bm(1,:)).^2 + (Jcons(2)*omegaB_Bm(2,:)).^2 + (1*omegaB_Bm(3,:)).^2);
    subplot(3,2,5); 
    plot(t,Ek/Ek0,'-','Color',col(:,iC)); hold on; grid on; ylabel('Normalized Variable');
    plot(t,h/h0,'--','Color',col(:,iC)); xlabel('Time [s]');
    
    % Modeled measurements
    rigidBodyRotation = inertiaRatiosOpt{i}.getRigidBodyRotation('opt');
%     RigidBodyRotation(Jopt(:,i),eye(3),omegaB_Bm(:,1),'omega0',0);
    omegaB_B = rigidBodyRotation.predictOmega(t-(i-1)*4*Tapprox);
    subplot(3,2,2); plot(t,omegaB_B(1,:),'-.k');  
    subplot(3,2,4); plot(t,omegaB_B(2,:),'-.k');  
    subplot(3,2,6); plot(t,omegaB_B(3,:),'-.k');     
    
end

subplot(3,2,2);
title('Body Angular Velocities'); legend('Unbiased, Aligned, Measurements','Local Fit');
subplot(3,2,5); 
legend('Normalized E_k','Normalized h'); title('Kinetic Energy and Angular Momentum');

    %}
