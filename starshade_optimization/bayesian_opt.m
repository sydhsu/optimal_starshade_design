%% Bayesian opt
clear all; close all; clc;

% Access origami pattern functions
addpath("crease_pattern_generator_package\");


% Seed and set random number generator
rng(0,'twister');

% Define variables and allowed lower/upper bounds
N = optimizableVariable('N',[6,30],'Type','integer');
n = optimizableVariable('n',[6,30],'Type','integer');
h = optimizableVariable('h',[1e-3,1e-2],'Type','real');
A = optimizableVariable('A',[1,5],'Type','real');
l = optimizableVariable('l',[1e-3,1e-1],'Type','real');
w = optimizableVariable('w',[1e-3,1e-1],'Type','real');
x = [N;n;h;A;l;w];

%x = [N,n,h,A,l,w];

% Initial point
A = 1;         % radius of inner polygon [m]
N = 4;         % number of sides of polygon
h = 1e-3;     % thickness of membrane [m]
n = 2;         % number of vertices along major fold
l = 1e-3;
w = 1e-3;
x0 = [N;n;h;A;l;w]; 

%% 
% Run the Bayesian optimization
results = bayesopt(@weightConstrained,x, ...
                   'Verbose',1, ...
                   'AcquisitionFunctionName','expected-improvement-plus', ...
                   'NumCoupledConstraints',11, ...
                   'PlotFcn', 'all',...
                   'MaxObjectiveEvaluations',60, ...
                   'IsObjectiveDeterministic',true,...
                   'ExplorationRatio', 0.8,...
                   'InitialX',array2table(x0','VariableNames',{'N','n','h','A','l','w'})); 

% Without constraints
% results = bayesopt(@weight,x, ...
%                    'Verbose',1, ...
%                    'AcquisitionFunctionName','expected-improvement-plus', ...
%                    'PlotFcn', 'all',...
%                    'MaxObjectiveEvaluations',60, ...
%                    'IsObjectiveDeterministic',true); 

% Verbose - how much output to print to the command line
% check constraints
% 'probability-of-improvement', 'lower-confidence-bound'
% 'ExplorationRatio', 0.5 default
% 'GPActiveSetSize',300 default


%% Best design obtained (single-objective opt.)
if ~isempty(results.XAtMinObjective)
    best = results.XAtMinObjective;
    
    N_b = best.N;
    n_b = best.n;
    h_b = best.h;
    A_b = best.A;
    l_b = best.l;
    w_b = best.w;
    
    % Visualization
    visualizeFlasher(N_b,n_b,h_b,A_b);
    best_weight = results.MinObjective;
    disp(['Weight of optimal design: ',best_weight])

else
    disp('No feasible design found!');
end


%% Helpers
function [wt,c] = weightConstrained(x)
    % Weight objective function that also returns constraints to be used
    % for bayesopt()
    % Input: x - design variables stored as table
    % Output: total weight of bars in structure

    rho_carbonfiber = 1550; % kg/m^3

    % Access from input table
    N = x.N;
    n = x.n;
    h = x.h;
    A = x.A;
    l = x.l;
    w = x.w;
    A_bar = l*w;
    
    % Construct the Flasher pattern
    [nodes_unfolded,~,edges,~,~,~] = flasher(N, n, h, A);

    % Bar lengths
    lengths = getEdgeLengths(nodes_unfolded, edges);
    
    % Objective function to minimize
    wt = sum(lengths*A_bar)*rho_carbonfiber;

    % Evaluate constraints to be returned as a vector
    x = [N;n;h;A;l;w];
    [c_ineq,c_eq] = constraintsByType(x);
    c = [c_ineq;c_eq];

end


