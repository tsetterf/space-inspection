function quats = integrateOmega( t, omega )
%INTEGRATEOMEGA Integrates angular velocity sequence and returns the quaternion
%   Based on Appendix D5 in [Farrell, 2008, Aided Navigation]
%   Inputs       
%       t       = 1xnT vector of times corresponding to the given omega_Bt vector [s]
%       omega   = 3xnT matrix of rotational velocities in the body frame [rad/s]
%   Outputs
%       quats   = 4xnT matrix of integrated unit quaternions, starting at [0 0 0 1]'

% Get length of data
nT = length(t);

% Allocate space for quaternions, and initialize the first one
quats = zeros(4,nT);
quats(:,1) = [0 0 0 1]';

% Integrate 
for i = 2:nT
   
    % Get time interval
    dt = t(i) - t(i-1);
    
    % Get a vector and matrix for integration
    % [Farrell2008 AppD5 with modification for our quaternion order and conventions]
    w = 1/2 * omega(:,i-1) * dt;
    W = [  -skew(w)     w       ;
           -w'          0       ];
    
    % Get the next quaternion [Farrell2008 AppD5]
    quats(:,i) = ( cos(norm(w))*eye(4) + sin(norm(w))/norm(w)*W ) * quats(:,i-1);
     
end


end

