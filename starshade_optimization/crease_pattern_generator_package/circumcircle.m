function [O, R] = circumcircle(r1, r2, r3)
    A = (r1 + r2)/2;
    B = (r2 + r3)/2;
    R90 = [0, 1; -1, 0];
    
    mat = [R90*(r2 - r1), -R90*(r3 - r2)];
    
    X = mat\(B - A);
    
    O1 = X(1)*R90*(r2 - r1) + A;
    O2 = X(2)*R90*(r3 - r2) + B;
    O = (O1 + O2)/2;
    
    R = (norm(O - r1) + norm(O - r2) + norm(O - r3))/3;   
end