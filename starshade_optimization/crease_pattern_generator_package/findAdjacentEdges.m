% Returns the four edges "adjacent" to the edge defined by edge0Index
% At node0, edge1 makes the smallest counterclockwise angle to edge0
%      and, edge2 makes the smallest clockwise angle to edge0
% At node1, edge3 makes the smallest counterclockwise angle to edge0
%      and, edge4 makes the smallest clockwise angle to edge0

function [edge1, edge2, edge3, edge4] = findAdjacentEdges(nodes, edges, edge0)    
    node0 = edges(edge0, 1);
    node1 = edges(edge0, 2);
    
    [edge1, edge2] = findAdjacentEdgesAtNode(nodes, edges, edge0, node0);
    [edge3, edge4] = findAdjacentEdgesAtNode(nodes, edges, edge0, node1);
end