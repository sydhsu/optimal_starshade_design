%% Flasher visualization

function visualizeFlasher(N,n,h,A)
    [nodes_unfolded, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasher(N,n,h,A);
    angles = foldedCreaseAngles(nodes_folded, edges, triangulated);
    % 2d folded plot
    plot2dNodesEdges(nodes_unfolded, edges, angles);
    %plotCreasePattern(nodes_unfolded,edges,.01)
    %plotCreasePatternWithLabels(nodes_unfolded, edges)

    % 3d folded plot
    plot3dNodesEdges(nodes_folded, edges, angles); rotate3d on; axis on; grid on;
%     f3 = plot3dNodesEdges(nodes_folded, edges, angles); rotate3d on;
%     plot3dCreasePatternWithLabels(nodes_folded, edges); rotate3d on; axis on;
%      
%     figure(f3);
%     p1 = patch();
%     p1.Faces = triangulated;
%     p1.Vertices = nodes_folded';
%     p1.FaceAlpha = 0.9;
%     p1.FaceColor = [0.7 0.7 0.7];
%     p1.LineStyle = 'none';
%     axis equal
%     axis off
%     view(3)
end
