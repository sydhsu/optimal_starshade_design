function [p_plus, p_minus] = intersectionOfThreeSpheres(p1, p2, p3, r1, r2, r3)
    d = norm(p2 - p1);
    
    i = dot(p3 - p1, p2 - p1)/norm(p2 - p1);
    j = norm(p3 - p1 - i*(p2 - p1)/norm(p2 - p1));
    
    ppx         = (r1^2 - r2^2 + d^2)/2/d;
    ppy         = (r1^2 - r3^2 + i^2 + j^2)/2/j - i*ppx/j;
    ppz_plus    = (r1^2 - ppx^2 - ppy^2)^(1/2);
    ppz_minus   = -ppz_plus;
    
    if not(isreal(ppz_plus))
        p_plus = nan(3, 1);
        p_minus = nan(3, 1);
        disp('cannot find 3 sphere intersection')
        return
    end
    
    pp_plus     = [ppx; ppy; ppz_plus];
    pp_minus    = [ppx; ppy; ppz_minus];     
    
    xp = (p2 - p1)/norm(p2 - p1);
    zp = cross(p2 - p1, p3 - p1)/norm(cross(p2 - p1, p3 - p1));
    yp = cross(zp, xp);
            
    p_plus      = [xp yp zp]*pp_plus + p1;
    p_minus     = [xp yp zp]*pp_minus + p1;
end