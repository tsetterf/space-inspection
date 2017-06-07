function [ q ] = rot2quat( R )
%ROT2QUAT Convert a rotation matrix to a quaternion in the SPHERES convention

% Note that Eigen in C++ will give you the conjugate of this

		div1 = 0.5*sqrt(1+R(1,1)+R(2,2)+R(3,3));
		div2 = 0.5*sqrt(1+R(1,1)-R(2,2)-R(3,3));
		div3 = 0.5*sqrt(1-R(1,1)-R(2,2)+R(3,3));
		div4 = 0.5*sqrt(1-R(1,1)+R(2,2)-R(3,3));
        
        q = [0 0 0 1];

        if (abs(div1) > 10e-5) 
			q(4) = div1;
			q(1) = 0.25*(R(2,3)-R(3,2))/q(4);
			q(2) = 0.25*(R(3,1)-R(1,3))/q(4);
			q(3) = 0.25*(R(1,2)-R(2,1))/q(4);
		elseif (abs(div2) > 10e-5) 
			q(1) = div2;
			q(2) = 0.25*(R(1,2)+R(2,1))/q(1);
			q(3) = 0.25*(R(1,3)+R(3,1))/q(1);
			q(4) = 0.25*(R(2,3)+R(3,2))/q(1);
        elseif ((abs(div3) > 10e-5))
			q(3) = div3;
			q(1) = 0.25*(R(1,3)+R(3,1))/q(3);
			q(2) = 0.25*(R(2,3)+R(3,2))/q(3);
			q(4) = 0.25*(R(1,2)-R(2,1))/q(3);
		else
			q(2) = div4;
			q(1) = 0.25*(R(1,2)+R(2,1))/q(2);
			q(3) = 0.25*(R(2,3)+R(3,2))/q(2);
			q(4) = 0.25*(R(3,1)-R(1,3))/q(2);
        end
        
        q = q' / norm(q);

end