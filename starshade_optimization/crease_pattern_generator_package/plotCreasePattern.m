% plots the 2D crease pattern, with a certain amount of radius left around
% the vertices to allow for better laser scoring

function fig = plotCreasePattern(nodes, edges, radius, fig)
    if nargin == 3
        fig = figure('Color', [1 1 1]);
    elseif nargin == 4
        figure(fig);
    end
        
    nEdges = size(edges,1);
    
    hold on
    
    for i = 1:nEdges        
        p1 = nodes(:,edges(i,1));
        p2 = nodes(:,edges(i,2));
        
        t = (p2-p1)/norm(p2-p1);
        p1prime = p1 + radius*t;
        p2prime = p2 - radius*t;
        
        line([p1prime(1) p2prime(1)],[p1prime(2) p2prime(2)],'Color',[0 0 0],'LineWidth',2)
    end
    
    %plot(nodes(1,:),nodes(2,:),'.')

    hold off
    
    %axis([-15 15 -15 15])
    axis equal
    %axis tight
    axis off
    
    
end