function [fig] = plotCreasePatternWithLabels(nodes, edges, fig)
    if nargin == 2    
        fig = plotCreasePattern(nodes, edges, 0);
    elseif nargin == 3
        plotCreasePattern(nodes, edges, 0, fig);
    end
    
    nNodes = size(nodes,2);
    figure(fig)
    hold on
    
    plot(nodes(1,:),nodes(2,:),'.','MarkerSize',3)
    
    labels = cellstr(num2str((1:nNodes)'));
    text(nodes(1,:),nodes(2,:),labels, 'VerticalAlignment', 'top',...
                                         'HorizontalAlignment', 'left',...                                        
                                         'BackgroundColor','w',...
                                         'FontSize',8,...
                                         'Margin',1,...
                                         'Color','k')
     
    
    hold off
end