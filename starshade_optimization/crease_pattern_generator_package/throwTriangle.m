function [p3_planar_plus, p3_planar_minus] = throwTriangle(p1_folded, p2_folded, p3_folded, p1_planar, p2_planar)
    b1 = (p2_folded - p1_folded)/norm(p2_folded - p1_folded); 
    a1 = (p2_planar - p1_planar)/norm(p2_planar - p1_planar);
    a2 = [0, -1, 0; 1, 0, 0; 0, 0, 1]*a1;
    
    r31p        = p3_folded - p1_folded;       
    r31p_dot_b1 = dot(r31p, b1);
    r31p_dot_b2 = norm(r31p - r31p_dot_b1*b1);
    
    r32p        = p3_folded - p2_folded;
    r32p_dot_b1 = dot(r32p, b1);
    r32p_dot_b2 = norm(r32p - r32p_dot_b1*b1);
    
    r31_plus        = r31p_dot_b1*a1 + r31p_dot_b2*a2;    
    r32_plus        = r32p_dot_b1*a1 + r32p_dot_b2*a2;
    p3_planar_plus  = mean([(p1_planar + r31_plus), (p2_planar + r32_plus)], 2);
    
    r31_minus       = r31p_dot_b1*a1 - r31p_dot_b2*a2;
    r32_minus       = r32p_dot_b1*a1 - r32p_dot_b2*a2;
    p3_planar_minus = mean([(p1_planar + r31_minus), (p2_planar + r32_minus)], 2);
end
