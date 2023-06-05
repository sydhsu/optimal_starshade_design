% Find all the nodes that connect to node0
function connectedNodeIndices = findConnectedNodes(edges, node0Index)    
    edgesFlag            = (edges == node0Index);
    edgesFlag            = edgesFlag(:,1) | edgesFlag(:,2);
    reducedEdges         = edges(edgesFlag,:);    
    connectedNodeIndices = reducedEdges(not(reducedEdges == node0Index));
end