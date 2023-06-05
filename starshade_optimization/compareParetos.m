%% Compare Pareto curve by overlaying plots
clear all; close all; clc

% Path to Pareto search data
load('data\pattern_search\ps_data.mat')
x_ps = x_pareto;
f_ps = f_pareto;

% Path to genetic algorithm data
load('data\ga\ga_data.mat')
x_ga = x_par;
f_ga = f_par;

%% Plot on the same graph
figure()
mkr_size = 10;
scatter(f_ps(:,1),f_ps(:,2),mkr_size); hold on;
scatter(f_ga(:,1),f_ga(:,2),mkr_size);
xlabel('Weight (kg)')
ylabel('(-) Deployed Area (m$^2$)')
title('\textbf{Pareto frontier}')
legend('Pattern search','Genetic algorithms','Interpreter','latex')
grid on

