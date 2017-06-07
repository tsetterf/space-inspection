function vec = vee(skewMatrix)
%SKEW Get a vectro from a skew symmetric cross product matrix

    vec = [ -skewMatrix(2,3) skewMatrix(1,3) -skewMatrix(1,2) ]';

end

