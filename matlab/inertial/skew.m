function skewMatrix = skew(vec)
%SKEW Get a skew symmetric cross product matrix from a vector

    skewMatrix = [  0          -vec(3)      vec(2);
                    vec(3)      0          -vec(1);
                   -vec(2)      vec(1)      0           ];

end

