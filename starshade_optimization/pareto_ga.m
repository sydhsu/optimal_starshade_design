%% Genetic algorithm multiobjective opt

clear all; close all; clc
set(0,'defaulttextInterpreter','latex')

% Access origami pattern functions
addpath("crease_pattern_generator_package\");

tic

%% Setup 
% Initial design (feasible) built by Arya et al
N = 14;         % number of sides of polygon
n = 14;         % number of vertices along major fold
h = 16e-3;      % thickness of membrane [m]
A = 0.734;      % radius of inner polygon [m]
l = 32e-3;      % length of cross-section
w = 1e-3;       % width of cross-section
x0 = [N,n,h,A,l,w];      % initial design

%%
nvars = length(x0); % dimension of x

% Objective function
multiobj_fn = @starshadeObjectives;

% Lower and upper bounds on design vars x = [N;n;h;A;l;w]
% Note all length units are in meters
lb = [3   2  1e-3 1e-3 1e-3 1e-3];
% ub = [30  30  0.1  5  0.1  0.1];
ub = [8  8  0.1  5  0.1  0.1];
A = []; b = []; Aeq = []; beq = [];
% Must be integers
intcon = [1 2];

options = optimoptions('gamultiobj', ...
                       'PlotFcn',{@gaplotpareto,@gaplotmaxconstr}, ... %@gaplotmaxconstr, @gaplotgenealogy, 'MaxGenerations',200*nvars, ...
                       'ConstraintTolerance',1e-4, ...
                       'FunctionTolerance',1e-4, ...
                       'PopulationSize',200, ...
                       'MaxTime',3600, ...
                       'MaxGenerations',100, ...
                       'InitialPopulationMatrix',x0 ... %'ParetoFraction', 0.35, ... % 0.35 default 
                       );

% Call the multiobjective optimization
[x_par,f_par,exitflag,output,population,scores] = gamultiobj(multiobj_fn,nvars,A,b,Aeq,beq,lb,ub,@constraintsByType,intcon,options);

runtime = toc; % end of computing the Pareto frontier
disp(['Elapsed time: ',num2str(runtime),' sec'])

%% Plot the Pareto frontier
figure()
scatter(f_par(:,1),f_par(:,2),5);
xlabel('Weight (kg)')
ylabel('(-) Deployed Area (m$^2$)')
title('\textbf{Pareto frontier}')
grid on

%% Design verification
x_pick = x_par(2,:);
% x_pick = x_ps(1,:);
%x_pick = [27	30	0.00254687500000000	3.67214062500000	0.00139022806586814	0.0409822228100349];

N = round(x_pick(1));
n = round(x_pick(2));
h = x_pick(3);
A = x_pick(4);
l = x_pick(5);
w = x_pick(6);

visualizeFlasher(N,n,h,A)

% figure()
% line([-w/2,-w/2,w/2,w/2,-w/2],[-l/2,l/2,l/2,-l/2,-l/2])
% axis equal

% Calculate Flasher specs
flasherSpec = analyzeFlasher(x_pick);
wt1 = flasherSpec(1);
area1 = flasherSpec(2);
R_stowed = flasherSpec(3);
h_stowed = flasherSpec(4);

fprintf('---------- Optimal design ---------- \n')
fprintf('Weight: %.4f \n',wt1)
fprintf('Deployed area: %.4f \n',area1)
fprintf('Stowed radius: %.4f \n',R_stowed)
fprintf('Stowed height: %.4f \n',h_stowed)

