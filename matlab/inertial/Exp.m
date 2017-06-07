function R = Exp( theta )
%EXP Creates a rotation matrix from an angle-axis vector using the SO(3) exponential map.
% Inputs:
%       theta   = 3x1 angle axis vector
% Outputs:
%       R       = 3x3 rotation matrix R = exp(theta^)

% Get skew-symmetric matrix
thetaCross = skew(theta);

% Get the magnitude of theta
thetaMag = norm(theta);

% Get rotation matrix
if thetaMag == 0
    R = eye(3);
else
    R = eye(3) + sin(thetaMag)/thetaMag * thetaCross + (1-cos(thetaMag))/thetaMag^2 * thetaCross^2;
end

end

