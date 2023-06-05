clear all; close all; clc;

% A   = 6e-2;
% N   = 9;
% h   = 5e-3;
% k   = 0.5;
% n   = 12;

% A   = 787.4/2;
% N   = 12;
% h   = 16;
% k   = 0.5;
% n   = 21;

% A   = 800;  % radius of inner polygon
% N   = 12;   % number of sides of polygon
% h   = 32;   % thickness
% n   = 11;   % number of minor folds

% A   = 400;
% N   = 8;
% h   = 50;
% k   = 0.5;
% n   = 7;

% A   = 400;  % radius of inner polygon
% N   = 6;   % number of sides of polygon
% h   = 10;   % thickness
% n   = 7;   % number of minor folds

A   = 2; % [m] %400;  % radius of inner polygon
N   = 6;   % number of sides of polygon
h   = 16e-3; % [m] thickness
n   = 14;

%[nodes_planar, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasher(N, n, h, A);
% conical flasher
phi_deg = 10; %5.45; % [deg]
% [nodes_planar, nodes_folded, edges, faces] = flasherConical_v3(N, n, h, A,phi_deg*pi/180);
[nodes_planar, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasherConical_v3(N, n, h, A,phi_deg*pi/180);

% Add last outer edge to make the pattern a closed shape (as-built Starshade OS test article)
% outerNodes = [2*n:n:2*N*n,n];
% outerEdges = reshape(outerNodes,[2,length(outerNodes)/2])';
% edges = [edges; outerEdges];
% outerPanels = [outerEdges,outerEdges(:,2)-1];
% triangulated = [triangulated;outerPanels];
% tri_faces = [tri_faces;outerPanels];

lengths     = getEdgeLengths(nodes_planar, edges);
lengths_p   = getEdgeLengths(nodes_folded, edges);
error       = (lengths - lengths_p)./lengths_p*100;

figure();
plot(1:length(edges), error, 'o--')
xlabel('Edge index')
ylabel('Length error (percent)')

angles = foldedCreaseAngles(nodes_folded, edges, triangulated);

% 2d folded plot
plot2dNodesEdges(nodes_planar, edges, angles);

% 3d folded plot
plot3dNodesEdges(nodes_planar, edges, angles); rotate3d on; axis on;
f3 = plot3dNodesEdges(nodes_folded, edges, angles); rotate3d on;
plotCreasePatternWithLabels(nodes_planar, edges)

plotCreasePattern(nodes_planar,edges,.01)

D_deployed  = 2*max(sqrt(nodes_planar(1, :).^2 + nodes_planar(2, :).^2));
% R_depl = D_deployed / 2
% R_deployed = max(vecnorm(nodes_planar,2,1));

%%
figure(f3);
p1 = patch();
p1.Faces = triangulated;
p1.Vertices = nodes_folded';
p1.FaceAlpha = 0.9;
p1.FaceColor = [0.7 0.7 0.7];
p1.LineStyle = 'none';
axis equal
axis off
view(3)

% figure();
% p2 = patch();
% p2.Faces = faces;
% p2.Vertices = nodes_planar';
% p2.FaceAlpha = 0.2;
% p1.FaceColor = [0.7 0.7 0.7];
% axis equal
% axis off 