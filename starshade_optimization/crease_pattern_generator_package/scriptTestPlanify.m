clear all; close all; clc;

R   = 20;
r   = 2;

za  = 2;
zb  = 0.5;

N   = 12;
n   = 12;

points3d = new_weird_gore(R, r, za, zb, N, n)

fig = figure();
plot3triangle(points3d, fig)
axis equal

points2d = planify(points3d)

fig = figure();
plot3triangle(points2d, fig)
axis equal
axis off

% n_points = 2*(n + 1);
% for i = 1:(n_points - 1)
%     norm(points3d(:, i + 1) - points3d(:, i)) - norm(points2d(:, i + 1) - points2d(:, i))
% end
