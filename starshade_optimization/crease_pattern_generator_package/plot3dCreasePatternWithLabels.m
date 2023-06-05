function [fig] = plot3dCreasePatternWithLabels(nodes, edges)    

    fig = plot3dNodesEdges(nodes, edges, nan(length(edges), 1));
    
    
    nNodes = size(nodes,2);
    figure(fig)
    hold on
    
    plot3(nodes(1,:), nodes(2,:), nodes(3,:), '.','MarkerSize',3)
    
    labels = cellstr(num2str((1:nNodes)'));
    text(nodes(1,:),nodes(2,:),nodes(3,:),labels, 'VerticalAlignment', 'top',...
                                         'HorizontalAlignment', 'left',...                                        
                                         'BackgroundColor','w',...
                                         'FontSize',8,...
                                         'Margin',1,...
                                         'Color','k')
     
    
    hold off
end