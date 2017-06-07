function theta = Log( R )
%LOG Gets an angle-axis vector from a rotation matrix using the SO(3) logarithmic map.
% Inputs:
%       R       = 3x3 rotation matrix
% Outputs:
%       theta   = 3x1 angle axis vector where -pi <= ||theta|| <= pi

% Get the magnitude of theta
thetaMag = acos( (trace(R)-1)/2 );

% Get the skew-symmetric matrix
thetaCross = thetaMag/(2*sin(thetaMag)) * (R - R');

% Get the angle-axis vector
theta = vee(thetaCross);

if sum(sum(R == eye(3))) == 9
    theta = [0 0 0]';
end

end

