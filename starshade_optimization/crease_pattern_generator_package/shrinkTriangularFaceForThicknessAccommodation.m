function nodes_shifted = shrinkTriangularFaceForThicknessAccommodation(nodes, folded_angles, thickness)    
    nodes_shifted = nan(2, 3);   
    
    node_perm = [   1 2 3;...
                    2 3 1;...
                    3 1 2];    
                
    angle_perm = [  1 3;...
                    2 1;...
                    3 2];    
    
    for i = 1:3        
        p   = nodes(:, node_perm(i, 1));
        p1  = nodes(:, node_perm(i, 2));
        p2  = nodes(:, node_perm(i, 3));
        
        alpha = abs(angleBetweenVectors(p2 - p, p1 - p));        
        
        rho1 = folded_angles(angle_perm(i, 1));
        rho2 = folded_angles(angle_perm(i, 2));                    
        
        if isnan(rho1)
            a1 = thickness*1e-6;            
        else            
            a1 = abs(thickness/2/tan(rho1/2));
            if a1 > thickness/2
                a1 = thickness/2;
            end
        end
        
        if isnan(rho2)
            a2 = thickness*1e-6;            
        else            
            a2 = abs(thickness/2/tan(rho2/2));
            if a2 > thickness/2
                a2 = thickness/2;
            end
        end
                
        nodes_shifted(:, i) = shrinkVertex(p, alpha, p1, p2, a1, a2);        
    end
end

% function nodes_shifted = shrinkTriangularFaceForThicknessAccommodation(face, nodes, edges, folded_angles, thickness)    
%     nodes_shifted = nodes;
%     
%     for node_id = face        
%         nodes_not_p     = face(not(face == node_id));
%         node_id_1       = nodes_not_p(1);
%         node_id_2       = nodes_not_p(2);
%         
%         p   = nodes(:, node_id);        
%         p1  = nodes(:, node_id_1);
%         p2  = nodes(:, node_id_2);
%         
%         alpha = abs(angleBetweenVectors(p2 - p, p1 - p));
%         
%         edge1 = findConnectingEdgeIndex(node_id, node_id_1, edges);
%         edge2 = findConnectingEdgeIndex(node_id, node_id_2, edges);
%         
%         rho1 = folded_angles(edge1);
%         rho2 = folded_angles(edge2);
%         
%         a1 = thickness/2/tan(rho1/2);        
%         a2 = thickness/2/tan(rho2/2);
%         
%         nodes_shifted(:, node_id) = shrinkVertex(p, alpha, p1, p2, a1, a2);
%     end
% end

function p_shifted = shrinkVertex(p, alpha, p1, p2, a1, a2)
    sin_alpha1 = sqrt((sin(alpha)^2)/((a2/a1 + cos(alpha))^2 + sin(alpha)^2));
    cos_alpha1 = sqrt(1 - sin_alpha1^2);
    
    e1 = (p1 - p)/norm(p1 - p);
    e2 = [0 -1; 1 0]*e1;
    
    if angleBetweenVectors(p2 - p, p1 - p) < 0   
        p_shifted = a1*cos_alpha1/sin_alpha1*e1 + a1*e2 + p;
    else
        p_shifted = a1*cos_alpha1/sin_alpha1*e1 - a1*e2 + p;
    end
end