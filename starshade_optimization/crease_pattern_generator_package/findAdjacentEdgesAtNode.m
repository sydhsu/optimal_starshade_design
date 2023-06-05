% Returns the two edges "adjacent" to the edge0 at node0
% At node0, edge1 makes the smallest counterclockwise angle to edge0
%      and, edge2 makes the smallest clockwise angle to edge0

function [edge1, edge2] = findAdjacentEdgesAtNode(nodes, edges, edge0, node0)

    edge0Nodes = edges(edge0, :); % the two nodes that edge0 connects    
    nodeOther  = edge0Nodes(not(edge0Nodes == node0)); % node that edge0 connects to node0
    
    if isempty(nodeOther)
        display('Umm... dude. edge0 does not hook up to node0.')
    end
   
    v0      = nodes(:, node0);
    vOther  = nodes(:, nodeOther);
    connectedNodeIndices = findConnectedNodes(edges, node0);
    connectedNodeIndices(connectedNodeIndices == nodeOther) = [];
    
    if isempty(connectedNodeIndices)
        edge1 = nan;
        edge2 = nan;
        return
    end
    
    if length(connectedNodeIndices) == 1 % if there's only one edge at node0 other than edge0
        edge1 = findConnectingEdgeIndex(node0, connectedNodeIndices, edges);
        %edge2 = nan;
        edge2 = edge1;
        return
    end
    
    % find all angles between edge0 and all other edges at v0
    theta = nan(1, length(connectedNodeIndices));
    for i = 1:length(connectedNodeIndices)
        iNode = connectedNodeIndices(i);
        theta(i) = angleBetweenVectors(nodes(:, iNode) - v0, vOther - v0);
    end    
       
    thetaPositive = theta;
    thetaPositive(theta < 0) = nan;
    
    thetaNegative = theta;
    thetaNegative(theta > 0) = nan;
    
    [~, minCCWThetaIndex] = min(thetaPositive);
    [~, minCWThetaIndex]  = max(thetaNegative);
        
    edge1OtherNodeIndex = connectedNodeIndices(minCCWThetaIndex);
    edge2OtherNodeIndex = connectedNodeIndices(minCWThetaIndex);
    
    edge1 = findConnectingEdgeIndex(node0, edge1OtherNodeIndex, edges);
    edge2 = findConnectingEdgeIndex(node0, edge2OtherNodeIndex, edges);
end