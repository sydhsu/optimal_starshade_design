function [fig] = plotCreasePatternWithFiducials(nodes,edges,vertexDrawRadius)
    fig = plotCreasePattern(nodes,edges,vertexDrawRadius);
    figure(fig)
    hold on
    
    %plot(nodes(1,:),nodes(2,:),'.','MarkerSize',3)
    
    % coordinates of the fiducial triangle    
    pFid = [1, 4, 4; 0, 1, -1]*1e-3/2;
        
    nNodes = size(nodes,2);
    for i = 1:nNodes
        % find all the nodes that connect to node i
        id = (edges == i);
        id = id(:,1) | id(:,2);
        connectingEdges = edges(id,:);
        connectingEdges(connectingEdges == i) = 0;
        connectingNodes = connectingEdges(:,1) + connectingEdges(:,2);
        
        % sort these nodes by counterclockwise angle
        nConnectingNodes = length(connectingNodes);        
        vectorsToConnectingNodes = nodes(:,connectingNodes) - nodes(:,i)*ones(1,nConnectingNodes);
        
        referenceVector = [1; 0; 0];
        
        % find the counterclockwise angles of the edge to the reference
        cosTheta = nan(nConnectingNodes,1);
        sinTheta = nan(nConnectingNodes,1);
        normal = [0; 0; 1];
        
        for j = 1:nConnectingNodes
            a = [vectorsToConnectingNodes(:,j); 0];
            cosTheta(j) = dot(referenceVector, a)/norm(a)/norm(referenceVector);
            sinTheta(j) = det([referenceVector'; a'; normal'])/norm(a)/norm(referenceVector);
        end              
        theta = atan2(sinTheta,cosTheta);
              
        % sort the connecting nodes by angle
        theta = sort(theta);                      
        psi = [theta(2:end) - theta(1:end-1); theta(1) - theta(end)];
        phi = theta + psi/2;
        phi(end) = phi(end) + pi();       
        
        for j = 1:nConnectingNodes
            C = [cos(phi(j)), -sin(phi(j)); sin(phi(j)), cos(phi(j))];            
            pFidRotated = C*pFid + nodes(:,i)*ones(1,3);
            fill(pFidRotated(1,:),pFidRotated(2,:),'k')
        end       
    end  
 
    
    hold off
end
