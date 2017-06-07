function [ qFixed ] = quatfix( q )
%quatfix Fix the quaternions, ensuring that they are continuous

n = size(q,2);
qFixed = zeros(4,n);
qFixed(:,1) = q(:,1);

for i=2:n
   qprev = q(:,i-1);
   qcurr = q(:,i);
   Rprev = quat2rot(qprev);
   Rcurr = quat2rot(qcurr);
   Rdiff = Rcurr * Rprev';
   qdiff = rot2quat(Rdiff);
   new_q = quatmult(qdiff,qprev);
   
   q(:,i) = new_q;
   qFixed(:,i) = new_q;
end


end

