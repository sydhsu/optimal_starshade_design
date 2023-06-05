function [p3_conical_plus, p3_conical_minus] = throwTriangleConical_v2(p1, p2, p3, p1_conical, p2_conical, phi)
    
    p = @(R, theta) [R*cos(theta); R*sin(theta); R*tan(phi)];
    fun = @(coords) [   norm(p(coords(1), coords(2)) - p1_conical) - norm(p3 - p1);...
                        norm(p(coords(1), coords(2)) - p2_conical) - norm(p3 - p2)];
    
    v_hat = (p2 - p1)/norm(p2 - p1);
    W = dot(p3 - p1, v_hat)*v_hat + p1;
    
    u_hat = (p2_conical - p1_conical)/norm(p2_conical - p1_conical);
    M = dot(p3 - p1, v_hat)*u_hat + p1_conical;
    
    theta_M = atan2(M(2), M(1));
    n_hat = [-tan(phi)*cos(theta_M); -tan(phi)*sin(theta_M); 1]*(1/sqrt(1 + tan(phi)^2));
    rho_hat_0 = cross(u_hat, n_hat);
    
    r_hat = [cos(theta_M); sin(theta_M); 0];
    
    if dot(rho_hat_0, r_hat) >= 0
        rho_hat = rho_hat_0;
    else
        rho_hat = -rho_hat_0;
    end
                    
    p3_conical_plus_0   = M + norm(p3 - W)*rho_hat;
    p3_conical_minus_0  = M - norm(p3 - W)*rho_hat;
    
    R_p3_conical_plus_0 = sqrt(p3_conical_plus_0(1)^2 + p3_conical_plus_0(2)^2);
    theta_p3_conical_plus_0 = atan2(p3_conical_plus_0(2), p3_conical_plus_0(1));
    
    R_p3_conical_minus_0 = sqrt(p3_conical_minus_0(1)^2 + p3_conical_minus_0(2)^2);
    theta_p3_conical_minus_0 = atan2(p3_conical_minus_0(2), p3_conical_minus_0(1));
    
    p3_conical_plus_coords     = fsolve(fun, [R_p3_conical_plus_0, theta_p3_conical_plus_0], optimoptions('fsolve', 'Display', 'off')); 
    p3_conical_minus_coords    = fsolve(fun, [R_p3_conical_minus_0, theta_p3_conical_minus_0], optimoptions('fsolve', 'Display', 'off'));
    
    p3_conical_plus = p(p3_conical_plus_coords(1), p3_conical_plus_coords(2));
    p3_conical_minus = p(p3_conical_minus_coords(1), p3_conical_minus_coords(2));

    err31_plus = norm(p3_conical_plus - p1_conical) - norm(p3 - p1);
    err32_plus = norm(p3_conical_plus - p2_conical) - norm(p3 - p2);
            
    err31_minus = norm(p3_conical_minus - p1_conical) - norm(p3 - p1);
    err32_minus = norm(p3_conical_minus - p2_conical) - norm(p3 - p2);
    
%     disp([  'v2 Err 3-1, +: ' num2str(err31_plus)...
%             ', Err 3-2, +: ' num2str(err32_plus)...
%             ', Err 3-1, -: ' num2str(err31_minus)...
%             ', Err 3-2, -: ' num2str(err32_minus)])
end