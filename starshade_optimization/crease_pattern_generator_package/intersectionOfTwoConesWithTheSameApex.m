% O         : apex
% a1        : axis of first cone
% theta1    : half-angle of first cone
% a2        : axis of second cone
% theta2    : half-angle of second cone
% L         : distance between p and O

function [p_plus, p_minus] = intersectionOfTwoConesWithTheSameApex(O, a1, theta1, a2, theta2, L)
    a1 = a1/norm(a1);
    a2 = a2/norm(a2);
    
    alpha       = acos(dot(a1, a2));
    
    vx          = cos(theta1);
    vy          = (cos(theta2) - cos(alpha)*vx)/sin(alpha);
    vz_plus     = (1 - vx^2 - vy^2)^(1/2);
    vz_minus    = -vz_plus;     
    
    pp_plus     = [vx; vy; vz_plus];
    pp_minus    = [vx; vy; vz_minus];
    
    if not(isreal(pp_plus))
        p_plus  = nan(3, 1);
        p_minus = nan(3, 1);
        disp('cannot find intersection of cones')
        return
    end
           
    xp = a1;
    zp = cross(a1, a2)/norm(cross(a1, a2));
    yp = cross(zp, xp);
    
    p_plus      = [xp yp zp]*(L*pp_plus) + O;
    p_minus     = [xp yp zp]*(L*pp_minus) + O;    
end