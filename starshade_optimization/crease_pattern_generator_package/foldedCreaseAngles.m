function angles = foldedCreaseAngles(nodes_folded, edges, faces)
    nEdges = size(edges, 1);    
    angles = nan(1, nEdges);
    
    for i = 1:nEdges
        p1_index = edges(i, 1);
        p2_index = edges(i, 2);
        
        faces_with_node1        = (faces(:, 1) == p1_index) | (faces(:, 2) == p1_index) | (faces(:, 3) == p1_index);
        faces_with_node2        = (faces(:, 1) == p2_index) | (faces(:, 2) == p2_index) | (faces(:, 3) == p2_index);
        faces_adjacent_to_edge  = faces((faces_with_node1 & faces_with_node2), :);
        
        if size(faces_adjacent_to_edge, 1) ~= 2
            %disp([num2str(size(faces_adjacent_to_edge, 1)) ' faces at edge connecting nodes ' num2str(p1_index) ' and ' num2str(p2_index)])            
            angles(i) = nan;
        else
            face1 = faces_adjacent_to_edge(1, :);
            face2 = faces_adjacent_to_edge(2, :);
            
            p3_1_index = face1((face1 ~= p1_index) & (face1 ~= p2_index));
            p3_2_index = face2((face2 ~= p1_index) & (face2 ~= p2_index));
                        
            p1 = nodes_folded(:, p1_index);
            p2 = nodes_folded(:, p2_index);
            p3_1 = nodes_folded(:, p3_1_index);
            p3_2 = nodes_folded(:, p3_2_index); 
            
            q1 = dot((p3_1 - p1), (p2 - p1))*(p2 - p1)/(norm(p2 - p1)^2) + p1;
            q2 = dot((p3_2 - p1), (p2 - p1))*(p2 - p1)/(norm(p2 - p1)^2) + p1;
            
            angles(i) = acos(dot((p3_1 - q1),(p3_2 - q2))/norm(p3_1 - q1)/norm(p3_2 - q2)); % this is the same as the angle between the two normals
           
            vec = mean([(p3_1 - q1) (p3_2 - q2)], 2);
            if vec(3) < 0 % this is a mountain
                angles(i) = -angles(i);
            end
        end
    end
end