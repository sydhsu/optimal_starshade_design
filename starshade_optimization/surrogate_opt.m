clc
clear all
close all

% Access origami pattern functions
addpath("crease_pattern_generator_package\");


% Seed and set random number generator
rng(0,'twister');

% Define variables and allowed lower/upper bounds

N = optimvar("N","LowerBound",6,"UpperBound",30,"Type","integer");
n = optimvar("n","LowerBound",6,"UpperBound",30,"Type","integer");
h = optimvar("h","LowerBound",1e-3,"UpperBound",1e-2);
A = optimvar("A","LowerBound",1,"UpperBound",5);
l = optimvar("l","LowerBound",1e-3,"UpperBound",.1);
w = optimvar("w","LowerBound",1e-3,"UpperBound",.1);
x = [N;n;h;A;l;w]';

fun = fcn2optimexpr(@getWeight,x,'OutputSize',[1,1]);
prob = optimproblem("Objective",fun);
opts = optimoptions("surrogateopt","MaxFunctionEvaluations",500);
rng(1,'twister') % For reproducibility
[sol,fval] = solve(prob,"Solver","surrogateopt","Options",opts)

% HELPERS
function [wt,c] = weightConstrained(x)
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
    [nodes_unfolded, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasher(N, n, h, A);

    % Bar lengths
    lengths = getEdgeLengths(nodes_unfolded, edges);
    
    % Objective function to minimize
    wt = sum(lengths*A_bar)*rho_carbonfiber;

    % Evaluate constraints to be returned as a vector
    x = [N;n;h;A;l;w];
    [c_ineq,c_eq] = constraintsByType(x);
    c = [c_ineq;c_eq];

end

function w = weight(x)
    % Input: x - design variables stored as table
    % Output: total weight of bars in structure

    rho_carbonfiber = 1550; % kg/m^3

    % Access from input table
    N = x.N;
    n = x.n;
    h = x.h;
    A = x.A;
    A_bar = x.l*x.w;
    
    % Construct the Flasher pattern
    % Prevent n from being 1
%     if n <= 1
%         w = Inf;
%     else
    [nodes_unfolded, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasher(N, n, h, A);
    
    % Deployed radius (in the x-y plane)
    % R_deployed  = max(vecnorm(nodes_unfolded,2,1)); 
    
    % Obtain the total surface area
    % TODO: function to compute area of triangle and quads
    
    % Obtain the number of bars in the structure
    % num_bars = size(edges,1);
    % Bar lengths
    lengths = getEdgeLengths(nodes_unfolded, edges);
    
    % Objective function to minimize
    w = sum(lengths*A_bar)*rho_carbonfiber;
%     end

end

function F = getWeight(x)
    rho_carbonfiber = 1550; % kg/m^3
    N = x(1);
    n = x(2);
    h = x(3);
    A = x(4);
    A_bar = x(5)*x(6);
    %x = [N;n;h;A;l;w]';
    % Construct the Flasher pattern
    if n <= 1
        w = Inf;
    else
    [nodes_unfolded, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasher(N, n, h, A);
    lengths = getEdgeLengths(nodes_unfolded, edges);
    w = sum(lengths*A_bar)*rho_carbonfiber;
    F.Fval = w;
    F.Ineq = constraintsByType(x);
    end

end