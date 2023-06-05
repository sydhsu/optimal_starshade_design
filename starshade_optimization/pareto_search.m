%% Pattern search multiobjective opt

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

nvars = length(x0); % dimension of x

% Objective function
multiobj_fn = @starshadeObjectives;

options = optimoptions('paretosearch', ...
                       'ParetoSetSize',200, ...
                       'PlotFcn',{@psplotparetof,@psplotfuncount,@psplotmaxconstr,}, ...
                       'InitialPoints',x0, ...
                       'ConstraintTolerance',1e-4, ...
                       'MaxTime',3600, ...
                       'ParetoSetChangeTolerance',1e-4, ...
                       'PollMethod','GPSPositiveBasis2N', ...
                       'MinPollFraction',1 ... % poll all points in the mesh
                       ); 

% Lower and upper bounds on design vars x = [N;n;h;A;l;w]
% Note all length units are in meters
lb = [3   2  1e-3 1e-3 1e-3 1e-3];
ub = [30  30  0.1  5  0.1  0.1];
A = []; b = []; Aeq = []; beq = [];

%% Execute Pareto search
[x_par,f_pareto,exitflag,output,residuals] = paretosearch(multiobj_fn,nvars,A,b,Aeq,beq,lb,ub,@constraintsByType,options);

% Round the integer design variables
x_pareto = x_par;
x_pareto(:,1) = round(x_pareto(:,1));
x_pareto(:,2) = round(x_pareto(:,2));

runtime = toc;
disp(['Elapsed time: ',num2str(runtime),' sec'])


%% Plot the Pareto frontier
figure()
scatter(f_pareto(:,1),f_pareto(:,2),'filled');
xlabel('Weight (kg)')
ylabel('(-) Deployed Area (m$^2$)')
title('\textbf{Pareto frontier}')
grid on


%% Design verification
x_pick = x_pareto(52,:);
% x_pick = [26	30	0.00103778366815476	3.76550099206349	0.00217032490079366	0.009];

N = round(x_pick(1));
n = round(x_pick(2));
h = x_pick(3);
A = x_pick(4);

visualizeFlasher(N,n,h,A)

% wt1 = weight(x_pick);
% area1 = deployedArea(x_pick);

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



%% Helpers

% function f = multiobj(x)
%     % Return the multiobjective function evaluation
%     f(:,1) = weight(x);
%     f(:,2) = -deployedArea(x); % want to minimize the negative = maximize
% end
% 
% function z = pickIndex(x,k)
%     z = multiobj(x); % evaluate both objectives
%     z = z(k); % return objective k (evaluation of kth function)
% end


% function w = weight(x)
%     % Unconstrained weight objective function
%     % Input: x - design variables stored as array
%     % Output: total weight of bars in structure
% 
%     rho_carbonfiber = 1550; % kg/m^3
% 
%     % Unpack variables
%     N = round(x(1));
%     n = round(x(2));
%     h = x(3);
%     A = x(4);
%     l = x(5); 
%     w = x(6);
%     A_bar = l*w;
%     
%     % Construct the Flasher pattern
%     [nodes_unfolded,~,edges,~,~,~] = flasher(N, n, h, A);
%     
%     lengths = getEdgeLengths(nodes_unfolded, edges);
%     
%     % Objective function to minimize
%     w = sum(lengths*A_bar)*rho_carbonfiber;
% 
% end
% 
% % Objective 2: deployed surface area
% function area = deployedArea(x)
%     % Input: x - design variables stored as array
%     % Output: total surface area in deployed state
% 
%     % Unpack variables
%     N = round(x(1));
%     n = round(x(2));
%     h = x(3);
%     A = x(4);
%     l = x(5); 
%     w = x(6);
% 
%     % Construct the Flasher pattern
%     [nodes_unfolded, ~, ~, ~, ~, ~] = flasher(N, n, h, A);
% 
% % Obtain the total surface area
% %     TODO: function to compute area of triangle and quads
% 
%     % Estimate area as a circumcircle with radius for now
%     R_deployed = max(vecnorm(nodes_unfolded,2,1));
%     area = pi*R_deployed^2;
% 
% 
% end

