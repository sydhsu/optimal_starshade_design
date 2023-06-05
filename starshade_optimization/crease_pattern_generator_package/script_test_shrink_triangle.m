clear all; close all; clc;

nodes = 5*rand(2, 3);
angles = pi()/4*rand(3, 1);
nodes_shifted = shrinkTriangularFaceForThicknessAccommodation(nodes, angles, 0.01);

fig = figure();

plot2triangle(nodes, fig)
plot2triangle(nodes_shifted, fig)

angles/pi()*180