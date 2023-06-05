function [p3_paraboloid_plus, p3_paraboloid_minus] = throwTriangleParaboloid(p1_folded, p2_folded, p3_folded, p1_conical, p2_conical, c)
    a_folded = p2_folded - p1_folded;
    a_folded = a_folded/norm(a_folded);
    
    q_folded = mean([   dot(p3_folded - p1_folded, a_folded)*a_folded + p1_folded,...
                        dot(p3_folded - p2_folded, a_folded)*a_folded + p2_folded], 2);
                    
    r = norm(p3_folded - q_folded);
            
    a_paraboloid = p2_conical - p1_conical;
    a_paraboloid = a_paraboloid/norm(a_paraboloid);
        
    q_paraboloid = mean([  dot(q_folded - p1_folded, a_folded)*a_paraboloid + p1_conical,...
                        dot(q_folded - p2_folded, a_folded)*a_paraboloid + p2_conical], 2);        
    
    fun = @(theta) norm(pHelper(theta, c, a_paraboloid, q_paraboloid) - q_paraboloid) - r;    
    
    theta = fsolve(fun, atan2(q_paraboloid(2), q_paraboloid(1)));    
    p3_paraboloid1 = pHelper(theta, c, a_paraboloid, q_paraboloid);
    
    p3_paraboloid1_antipodal = 2*q_paraboloid - p3_paraboloid1;        
    theta = fsolve(fun, atan2(p3_paraboloid1_antipodal(2), p3_paraboloid1_antipodal(1)));
    p3_paraboloid2 = pHelper(theta, c, a_paraboloid, q_paraboloid);
    
    % project onto a plane
    a1_projected = p2_conical(1:2) - p1_conical(1:2);
    a1_projected = a1_projected/norm(a1_projected);    
    
    a2_projected = [0, -1; 1, 0]*a1_projected;
    
    if dot(p3_paraboloid1(1:2) - q_paraboloid(1:2), a2_projected) > 0
        p3_paraboloid_plus = p3_paraboloid1;
        p3_paraboloid_minus = p3_paraboloid2;
    else
        p3_paraboloid_plus = p3_paraboloid2;
        p3_paraboloid_minus = p3_paraboloid1;        
    end
end

function p = pHelper(theta, c, a, q)
%     c   = 1/tan(phi);
%     p_z = dot(a, q)/dot(a, [c*cos(theta); c*sin(theta); 1]);    
%     p   = p_z*[c*cos(theta); c*sin(theta); 1];   

    r_cand = roots([a(3)*c, (a(1)*cos(theta) + a(2)*sin(theta)), -dot(a, q)]);
    
    %r = r_cand(isreal(r_cand))
    r = r_cand(end);
    
    %r = dot(a, q)/dot([cos(theta); sin(theta); tan(phi)], a);
    p = [r*cos(theta); r*sin(theta); c*r^2];
    
%     fun = @(r) dot([r*cos(theta); r*sin(theta); r*tan(phi)], a) - dot(q, a);
%     r = fsolve(fun, 0);
%     p = r*[cos(theta); sin(theta); tan(phi)];
end