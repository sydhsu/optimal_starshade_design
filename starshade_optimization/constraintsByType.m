function [c, ceq] = constraintsByType(x)
    % Returns a vector of the inequality constraints and a vector of the
    % equality constraints
    % Input: x - design vector [N;n;h;A;l;w]
    % Outputs: c_ineq - vector of evaluated inequality constraints
    %          c_eq - vector of evaluated equality constraints

    % Constant and uniform force
    F_max = 1e5; % [N] WAG; the effect is just to distribute this across all bars in the structure

    % Material properties of T800S carbon fiber composite
    E = 163e9; % Young's modulus [Pa]
    sigma_t = 3290e6; % material tensile yield strength [Pa]
    sigma_c = 1490e6; % material compressive yield strength [Pa]

    % Unpack the design variables
    N = round(x(1));
    n = round(x(2));
    h = x(3);
    A = x(4);
    l = x(5);
    w = x(6);

    % Cross-sectional bar area (rectangular)
    A_bar = l*w;

    % Construct the Flasher pattern
    if n <=1 
        c = ones(6,1)*Inf; % infinite penalty;
        ceq = [];
    else
    [nodes_unfolded, nodes_folded, edges, ~, ~, ~] = flasher(N, n, h, A);

    %%%%%%%%%
    %Estimate area as a circumcircle with radius for now
    R_deployed = max(vecnorm(nodes_unfolded,2,1));
    deployedArea = pi*R_deployed^2;
    %%%%%%%%%

    % Bar lengths
    lengths = getEdgeLengths(nodes_unfolded, edges);

    % Assume that the applied force is distributed evenly across all bars
    % as a first cut at analysis (complex to model in reality)
    F_app = F_max / (N*n);

    % Stowed radius of outermost panel (in the x-y plane)
    max_R_stowed = 5; % [m] WAG
    R_stowed = max(vecnorm(nodes_folded,2,1));

    % Axial (tensile and compressive) stress
    axial_stress_max = abs(F_app)/A_bar;
    FS_t = sigma_t / axial_stress_max;

    FS_c = sigma_c / axial_stress_max;

    % Maximum bending stress
    % First determine which dimension is smaller to determine the smallest
    % moment of inertia, and therefore the highest bending stress
    b = w; % set the smaller dimension
    h = l; % set the larger dimension
    if w>l
        b = l; % smaller
        h = w; % larger
    end

    I_min = 1/12*b^3*h; % smallest moment of inertia
%     y_max_bend = b/2; % location of maximum bending at the surface
%     bending_stress_max = Mmax*y_max_bend/I_min;
%     FS_bend = sigma_c / bending_stress_max;

    % Buckling stress
    % Calculate for every beam length in the structure
    buckling_stress_cr = max(pi^2*E*I_min ./ (lengths.^2) / A_bar); % pi^2*E*b^2/(12*lengths.^2)
%     y_max_bkl = b/2; % location of maximum buckling load
    %buckling_stress_max = Mmax*y_max_bkl/I_min; % applied buckling stress
    buckling_stress_max = abs(F_app/A_bar);
    FS_bkl = buckling_stress_cr / buckling_stress_max; % we want the critical load to be higher than the applied laod

    % Compile into constraint vector
    % Format: c <= 0
    mm = 1e-3; % 1 millimeter

%     c = [R_stowed - max_R_stowed;
%          -deployedArea+1000;
%          -FS_t + 1.4;
%          -FS_c + 1.4;
%          -FS_bkl + 1.4; 
%          -N + 3;
%          -n + 2;
%          -h + mm;
%          -A + mm;
%          -l + mm;
%          -w + mm;
%          ];

    c = [R_stowed - max_R_stowed;
         -deployedArea;
         -FS_t + 1.4;
         -FS_c + 1.4;
         -FS_bkl + 1.4; 
         -N + 3;
         -n + 2;
         -h + mm;
         -A + mm;
         -l + mm;
         -w + mm;
         ];

    ceq = []; 
    end


end