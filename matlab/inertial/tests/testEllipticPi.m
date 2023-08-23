%% Trapezoidal integration of Hurtado's integral as well as using ellipticPi to see equivalence

addpath('../');

% Test parameters
gammaBar = [1.239 1.1905]';
RB0toW = eye(3);
omegaB0_W = [0.699809513728002 0.716087962216080 0.355816171813249]';

% Create a test rigid body
ratios = InertiaRatios(gammaBar);
rigidBodyRotation = RigidBodyRotation(ratios,omegaB0_W,RB0toW);

% Times
t = 0:0.05:4*rigidBodyRotation.T;
nT = length(t);

% Setup variable histories
RBtoW = zeros(3,3,nT);          RBtoW(:,:,1) = RB0toW;
omegaB_B = zeros(3,nT);         omegaB_B(:,1) = RB0toW' * omegaB0_W;

% Get variables for ellipticPi integration
J = rigidBodyRotation.inertiaRatios.matrix();
J1 = J(1,1); J2 = J(2,2); J3 = J(3,3);
omega1m = rigidBodyRotation.omegaMax(1); 
omega2m = rigidBodyRotation.omegaMax(2); 
omega3m = rigidBodyRotation.omegaMax(3);
n = (J2^2*omega2m^2 - J3^2*omega3m^2) / (J3^2*omega3m^2);

%% Calculate Hurtado's alpha elliptic integral numerically
[sn,cn,dn] = ellipj(rigidBodyRotation.omegaP*t, rigidBodyRotation.k^2);
ellintNum = zeros(1,nT);
for i = 1:nT
    ellintNum(i) = trapz(t(1:i),1./(1+n*sn(1:i).^2));
end

% Calculate the integrand for plotting
integrandNum = 1./(1+n*sn.^2);

%% Calculate Hurtado's alpha elliptic integral analytically

% Calculate the quarter-period of the alpha elliptic integral in terms of u
Tu = rigidBodyRotation.T * rigidBodyRotation.omegaP;

% Determine number of revolutions for each time
numRev = floor(rigidBodyRotation.omegaP*t./(4*Tu));

% Get elliptic functions
[snUm,cnUm,~] = ellipj(rigidBodyRotation.omegaP*t,rigidBodyRotation.k^2);

% Use quadrants to properly determine phiM. Note that asin returns a value between [-pi/2 and pi/2]
phiM = zeros(1,length(t));
quadrant = zeros(1,length(t));
for i = 1:length(t)
    if snUm(i) >= 0 && cnUm(i) >= 0     % Q1
        phiM(i) = numRev(i)*2*pi + asin(snUm(i));
        quadrant(i) = 1;
    elseif snUm(i) >= 0 && cnUm(i) < 0  % Q2
        phiM(i) = numRev(i)*2*pi + pi - asin(snUm(i)); 
        quadrant(i) = 2;
    elseif snUm(i) < 0 && cnUm(i) < 0  % Q3
        phiM(i) = numRev(i)*2*pi + pi - asin(snUm(i));
        quadrant(i) = 3;
    else                               % Q4
        phiM(i) = numRev(i)*2*pi + 2*pi + asin(snUm(i));
        quadrant(i) = 4;
    end
end

% phiM = asin(snUm);
a0 = (J2^2*J3 - J2*J3^2)/(J1*J2 - J1*J3); 
ellintAna = 1/rigidBodyRotation.omegaP*ellipticPi(-n,phiM,rigidBodyRotation.k^2);

% Calculate the integrand used by the ellipticPi function
% theta = 0:0.01:phiM(end);
% integrandAna = 1./((1+n*sin(theta).^2).*sqrt(1-rigidBodyRotation.k^2*sin(theta).^2));

% Get the difference between the analytic and numerical methods of elliptic integration
ellintDiff = ellintAna - ellintNum;

% Plot
figure(1); clf;
subplot(4,1,1);
plot(t,ellintNum,'b'); hold on; plot(t,ellintAna,'--r'); grid on; title('Integral \int_0^t 1/(1+n*sn^2)');
legend('Numerical','Analytic');
subplot(4,1,2);
plot(t,integrandNum,'.-b'); title('Numerical Integrand 1/(1+n*sn^2)'); grid on;
subplot(4,1,3);
%plot(theta,integrandAna,'.-r'); title('Analytic Integrand 1./((1+n*sin(theta)^2)*sqrt(1-k^2*sin(theta)^2))');
plot(t,ellintDiff,'-r'); title('Difference between analytic and numerical integrals'); grid on; hold on;
subplot(4,1,4);
plot(t,quadrant,'.-r'); title('Quadrant'); grid on;
