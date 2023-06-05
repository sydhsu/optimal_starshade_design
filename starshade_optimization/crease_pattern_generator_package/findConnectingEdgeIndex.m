% finds the index of the edge that connects the nodes defined by node1 and
% node2

function edgeIndex = findConnectingEdgeIndex(node1, node2, edges)

    edgeIndices = 1:size(edges,1);
    
    node1Flag = (edges == node1);
    node1Flag = node1Flag(:,1) | node1Flag(:,2);
    node2Flag = (edges == node2);
    node2Flag = node2Flag(:,1) | node2Flag(:,2);
    
    edgeIndex = edgeIndices(node1Flag & node2Flag);

end