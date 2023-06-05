% plots in 3d the origami represented by nodes, edges

function fig = plot3dNodesEdges(nodes, edges, angles)
    fig = figure('Color',[1 1 1]);
    
    nEdges = size(edges,1);
        
    hold on
    for i = 1:nEdges        
        if isnan(angles(i))
            color = [0 0 0];        
        elseif angles(i) > 0
            color = [1 0 0] + angles(i)/pi()*[0 1 1];
        else
            color = [0 0 1] - angles(i)/pi()*[1 1 0];
        end
        
        %color = 'k';
        
        p1 = nodes(:,edges(i,1));
        p2 = nodes(:,edges(i,2));
        line([p1(1) p2(1)],[p1(2) p2(2)],[p1(3) p2(3)], 'Color', color, 'LineWidth', 1)
    end
    
    plot3(nodes(1,:), nodes(2,:), nodes(3,:), 'LineStyle', 'none', 'Marker', '.', 'MarkerEdgeColor', 'k')
       
    hold off
    
    view(38,30)
    axis equal
    axis vis3d
    axis off
    
end