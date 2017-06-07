function qnew = quatmult( q1, q2 )
%QUATMULT Multiply quaternions consistent with SPHERES convention and rotation matrix order
%   Inputs
%       q1      a column vector quaternion
%       q2      a column vector quaternion or matrix of several column vector quaternions
%   Outputs
%       qnew    the result of quaternion composition between q1 and column vectors of q2
%   Note     
%       The order is the same as rotation matrices 
%       (i.e. RCtoA = RBtoA * RCtoB <=> qCtoA = quatmult(qBtoA,qCtoB)

qnew = quatbarmat(q1) * q2;

end