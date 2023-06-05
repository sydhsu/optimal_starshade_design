% Draws a polygon defined by the vertices v_i = x_i,y_i, with the vertices
% rounded with radius

function roundedPolygon(x, y, radius, fig)
    n_vertices = length(x);
    if length(y) ~= n_vertices
        display('length(x) not equal to length(y)')
    end    

    line_x = nan;
    line_y = nan;

    for i = 1:n_vertices        
        if i == 1
            v_im1 = [x(end) y(end)];
        else        
            v_im1 = [x(i - 1) y(i - 1)];
        end
        
        v_i = [x(i) y(i)];
        
        if i == n_vertices
            v_ip1 = [x(1) y(1)];
        else        
            v_ip1 = [x(i + 1) y(i + 1)];
        end
        
        theta = angleBetweenVectors(v_ip1 - v_i, v_im1 - v_i);
        l     = radius/abs(tan(theta/2));
        d     = sqrt(radius^2 + l^2);

        pMid1 = (v_im1 + v_i)/2;
        pMid2 = (v_ip1 + v_i)/2;

        pEnd1 = (v_im1 - v_i)/norm(v_im1 - v_i)*l + v_i;
        pEnd2 = (v_ip1 - v_i)/norm(v_ip1 - v_i)*l + v_i;

        center_vec_dir = (v_im1 - v_i)/norm(v_im1 - v_i) + (v_ip1 - v_i)/norm(v_ip1 - v_i);
        center_vec     = center_vec_dir/norm(center_vec_dir)*d + v_i;    

        v_im1 = pEnd1 - center_vec;
        v_ip1 = pEnd2 - center_vec;

        theta1      = atan2(v_im1(2), v_im1(1));
        theta2      = theta1 + angleBetweenVectors(v_im1, v_ip1);
        d_theta     = 2*pi()/360*5;
        n_circle    = ceil(abs(theta2 - theta1)/d_theta);
        [Cx, Cy]    = create_open_circle(radius, center_vec(1), center_vec(2), theta1, theta2, n_circle);
        
        line_x      = [line_x pMid1(1) pEnd1(1) Cx pEnd2(1) pMid2(1)];
        line_y      = [line_y pMid1(2) pEnd1(2) Cy pEnd2(2) pMid2(2)];
    end
    
    figure(fig)
    line(line_x, line_y, 'Color', 'k')
end