function f = starshadeObjectives(x)
    % Objective functions to minimize
    % Input: x - design variables stored as array
    % Outputs: f - row vector of objective function evaluations: 
    % Weight, -Area

    rho_carbonfiber = 1550; % kg/m^3

    % Unpack variables
    N = round(x(1)); % must be an integer
    n = round(x(2));
    h = x(3);
    A = x(4);
    l = x(5); 
    w = x(6);
    A_bar = l*w;
    
    % Construct the Flasher pattern
    [nodes_unfolded,~,edges,~,~,~] = flasher(N, n, h, A);
    
    lengths = getEdgeLengths(nodes_unfolded, edges);
    
    % 1. Total weight of all bars (neglect weight of blanket)
    weight = sum(lengths*A_bar)*rho_carbonfiber;

    % 2. Estimate area as a circumcircle with radius for now
    R_deployed = max(vecnorm(nodes_unfolded,2,1));
    deployedArea = pi*(R_deployed^2);


    % Return as a single vector
    f(:,1) = weight;
    f(:,2) = -deployedArea; % negative because we want to maximize

end