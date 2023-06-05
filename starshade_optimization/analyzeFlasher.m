function f = analyzeFlasher(x)
    % Consolidate Flasher geometric analysis into a single multiobjective function
    % Input: x - design variables stored as array
    % Outputs: weight - total weight of bars in structure
    %          deployedArea - total surface area of deployed structure
    %          R_stowed - radius of stowed configuration
    %          h_stowed - height of stowed configuration

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
    [nodes_unfolded,nodes_folded,edges,~,~,~] = flasher(N, n, h, A);
    
    lengths = getEdgeLengths(nodes_unfolded, edges);
    
    % 1. Total weight of all bars (neglect weight of blanket)
    weight = sum(lengths*A_bar)*rho_carbonfiber;

    % 2. Estimate area as a circumcircle with radius for now
    R_deployed = max(vecnorm(nodes_unfolded,2,1));
    deployedArea = pi*R_deployed^2;

    % 3. Stowed radius
    R_stowed = max(vecnorm(nodes_folded,2,1));

    % 4. Stowed height
    h_stowed = max(nodes_folded(3,:));


    % Return as a single vector
    f(:,1) = weight;
    f(:,2) = deployedArea; % negative because we want to maximize
    f(:,3) = R_stowed;
    f(:,4) = h_stowed;

end