close all; clear all; clc;

for i = 1:100

    p1 = rand(3, 1);
    p2 = rand(3, 1);
    p3 = rand(3, 1);

    displacement = 0.01*rand(3, 1);
    %displacement = zeros(3, 1);

    p1p = p1 + displacement;
    p2p = p2 + displacement;

    z0 = rand(1);

    [p3p_plus, p3p_minus] = positionTriangleByZ(p1, p2, p3, p1p, p2p, z0);

end

% plot3dNodesEdges([p1, p2, p3, p1p, p2p, p3p_plus, p3p_minus], [1 2; 2 3; 3 1; 4 5; 5 6; 6 4; 5 7; 7 4], nan(1, 8))
% axis on
% grid on
% 
% norm(p3p_plus - p1p) - norm(p3 - p1)
% norm(p3p_plus - p2p) - norm(p3 - p2)
% 
% norm(p3p_minus - p1p) - norm(p3 - p1)
% norm(p3p_minus - p2p) - norm(p3 - p2)