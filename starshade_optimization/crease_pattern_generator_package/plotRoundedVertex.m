% plots a single rounded vertex of a polygon
% being part of the polygon, this vertex can only have degree 2
% so it also draws one half of each of the two edges connected to it
function plotRoundedVertex(v0, v1, v2, radius, fig)
    theta = angleBetweenVectors(v2 - v0, v1 - v0);
    l     = radius/abs(tan(theta/2));
    d     = sqrt(radius^2 + l^2);
    
    pMid1 = (v1 + v0)/2;
    pMid2 = (v2 + v0)/2;
    
    pEnd1 = (v1 - v0)/norm(v1 - v0)*l + v0;
    pEnd2 = (v2 - v0)/norm(v2 - v0)*l + v0;
          
    center_vec_dir = (v1 - v0)/norm(v1 - v0) + (v2 - v0)/norm(v2 - v0);
    center_vec     = center_vec_dir/norm(center_vec_dir)*d + v0;    
        
    v1 = pEnd1 - center_vec;
    v2 = pEnd2 - center_vec;
    
    theta1 = atan2(v1(2), v1(1));
    theta2 = theta1 + angleBetweenVectors(v1, v2);
       
    d_theta     = 2*pi()/360*5;
    n_circle    = ceil(abs(theta2 - theta1)/d_theta);
    [Cx, Cy]    = create_open_circle(radius, center_vec(1), center_vec(2), theta1, theta2, n_circle);
    
    figure(fig)
    line([pMid1(1) pEnd1(1) Cx pEnd2(1) pMid2(1)], [pMid1(2) pEnd1(2) Cy pEnd2(2) pMid2(2)], 'Color', 'k')    
end