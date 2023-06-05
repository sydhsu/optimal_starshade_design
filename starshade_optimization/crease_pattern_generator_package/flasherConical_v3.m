% Outputs the crease pattern (both unfolded and folded) for a
% thickness-accommodating flasher with N gores (or 2N major folds), n nodes
% per major fold line, h is the thickness, A is the radius of the hub
% ======================================================================= %
% Inputs:
%   N - number of sides of polygon
%   n - number of nodes per major fold line
%   h - panel thickness
%   A - radius of polygonal hub
%
% Outputs:
%   nodes_conical - 3-D coordinates in unfolded crease pattern
%   nodes_folded - 3-D coordinates in folded crease pattern
%   edges - array of node pairs that form an edge
%   triangulated - groups of 3 nodes after triangulating the entire pattern
%   tri_faces - groups of 3 nodes that form a triangular panel
%   quad_faces - groups of 4 nodes that form a quadrilateral panel
% ======================================================================= %

function [nodes_conical, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasherConical_v3(N, n, h, A, phi)
    beta    = 2*pi()/N;
    delta_z = 2*A*sin(beta/2)*tan(beta/2); % changes in height of mountain for a zero-thickness flasher
    d       = h/cos(beta/2);
    rot     = [ cos(beta), -sin(beta), 0;...
                sin(beta), cos(beta), 0;...
                0, 0, 1];

    m_conical   = nan(3, n);
    m_folded    = nan(3, n);
    v_conical   = nan(3, n);
    v_folded    = nan(3, n);

    m_conical(:, 1)     = [A; 0; A*tan(phi)];
    m_folded(:, 1)      = [A; 0; 0];
    v_conical(:, 1)     = m_conical(:, 1);
    v_folded(:, 1)      = m_folded(:, 1);
    v_folded(:, 2)      = rot*[A + 2*d; 0; 0];

    [m_conical(:, 2), m_folded(:, 2), v_conical(:, 2)] = ...
        flasherConicalSolveFirstPoint(m_conical, m_folded, v_conical, v_folded, rot, A, d, delta_z, phi);
    for i = 2:(n - 1)
        [m_conical(:, i + 1), m_folded(:, i + 1), v_conical(:, i + 1), v_folded(:, i + 1)] = ...
            flasherConicalSolveGeneralPoint(m_conical, m_folded, v_conical, v_folded, i, rot, A, d, delta_z, phi);
    end
    
    nodes_conical   = [];
    nodes_folded    = [];
    edges           = [];
    triangulated    = [];
    tri_faces       = [];
    quad_faces      = [];

    for j = 1:N
        nodes_conical   = [nodes_conical, (rot^(j - 1))*v_conical, (rot^(j - 1))*m_conical];
        nodes_folded    = [nodes_folded, (rot^(j - 1))*v_folded, (rot^(j - 1))*m_folded];
    end
    
    i_grid = (ones(N, 1)*(1:n))';   % node id along the gore
    j_grid = ones(n, 1)*(1:N);      % node id across the circle
    
    v_id = n*(2*j_grid - 2) + i_grid;
    m_id = n*(2*j_grid - 1) + i_grid;  
    v_id = [v_id, v_id(:, 1)];
    m_id = [m_id, m_id(:, 1)];   
   
    for j = 1:N 
        edges = [edges;...                
                v_id(1, j),     v_id(1, j + 1);...  % side of the polygon
                v_id(1, j),     v_id(2, j);...      % major valley
                v_id(1, j),     m_id(2, j);...      % major mountain
                m_id(2, j),     v_id(1, j + 1);...  % minor valley
                m_id(2, j),     v_id(2, j);...      % minor mountain                    
%                 m_id(2, j),     v_id(2, j + 1)];    % diagonal
        ];
                            
        triangulated = [triangulated;...
                v_id(1, j),     m_id(2, j), v_id(2, j);...                    
                v_id(1, j),     m_id(2, j), v_id(1, j + 1);...
                v_id(1, j + 1), m_id(2, j), v_id(2, j + 1)];

        tri_faces = [tri_faces;...
                v_id(1, j),     m_id(2, j), v_id(2, j);...                    
                v_id(1, j),     m_id(2, j), v_id(1, j + 1)];
                
        for i = 2:(n - 1)            
            edges = [edges;...
                    v_id(i, j),     v_id(i + 1, j);...  % major valley
                    m_id(i, j),     m_id(i + 1, j);...  % major mountain
                    m_id(i + 1, j), v_id(i, j + 1);...  % minor valley
                    m_id(i + 1, j), v_id(i + 1, j);...  % minor mountain
                    %m_id(i, j),     v_id(i + 1, j);...  % diagonal     
%                     v_id(i, j),     m_id(i + 1, j);... % cross-diagonal to one above            
%                     m_id(i + 1, j), v_id(i + 1, j + 1)];% diagonal   
            ];
                
            triangulated = [triangulated;...
                    %v_id(i, j),     m_id(i, j),     v_id(i + 1, j);...
                    %m_id(i, j),     v_id(i + 1, j), m_id(i + 1, j);...
                    v_id(i, j),     m_id(i + 1, j), m_id(i, j);...
                    v_id(i, j),     v_id(i + 1, j), m_id(i + 1, j);...
                    m_id(i, j),     m_id(i + 1, j), v_id(i, j + 1);...
                    v_id(i, j + 1), v_id(i + 1, j + 1), m_id(i + 1, j)];
            
            quad_faces = [quad_faces;...
                    v_id(i, j),     m_id(i, j),     m_id(i + 1, j),  v_id(i + 1, j);... 
                    m_id(i, j),     v_id(i-1,j+1),  v_id(i,j+1),     m_id(i+1, j)];
        end   
        edges = [edges;[m_id(n,j), v_id(n,j+1)]]; % adds outermost edge to close the outer polygon shape
    end
end

function [m_conical_2, m_folded_2, v_conical_2] = flasherConicalSolveFirstPoint(m_conical, m_folded, v_conical, v_folded, rot, A, d, delta_z, phi)
    fun = @(z2) helperFirst(m_conical, m_folded, v_conical, v_folded, rot, A, d, z2, phi);    
    mz2 = fsolve(fun, delta_z, optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12));    
    [m_conical_2, m_folded_2, v_conical_2] = ...
        generateCandidatesFirst(m_conical, m_folded, v_conical, v_folded, rot, A, d, mz2, phi);   
end
function [m_conical_2, m_folded_2, v_conical_2] = generateCandidatesFirst(m_conical, m_folded, v_conical, v_folded, rot, A, d, mz2, phi)    
    m_folded_2 = rot*[A + d; 0; mz2];    
    [m_conical_2, ~] = throwTriangleConical_v2(v_folded(:, 1), rot*v_folded(:, 1), m_folded_2, v_conical(:, 1), rot*v_conical(:, 1), phi);    
    [v_conical_2, ~] = throwTriangleConical_v2(m_folded(:, 1), m_folded_2, v_folded(:, 2), m_conical(:, 1), m_conical_2, phi);
end
function diff = helperFirst(m_conical, m_folded, v_conical, v_folded, rot, A, d, mz2, phi)
    [m_conical_2, m_folded_2, v_conical_2] = generateCandidatesFirst(m_conical, m_folded, v_conical, v_folded, rot, A, d, mz2, phi);    
    diff = norm(m_conical_2 - rot*v_conical_2) - norm(m_folded_2 - rot*v_folded(:, 2));    
end

function [m_conical_n, m_folded_n, v_conical_n, v_folded_n] = flasherConicalSolveGeneralPoint(m_conical, m_folded, v_conical, v_folded, i, R, A, d, delta_z, phi) 
    fun = @(z_n) helperGeneral(m_conical, m_folded, v_conical, v_folded, i, R, A, d, z_n(1), z_n(2), phi);
    mz_prev = m_folded(3, i);
    vz_prev = v_folded(3, i);
    z_n = fsolve(fun, [mz_prev + delta_z, vz_prev], optimoptions('fsolve', 'Display', 'off', 'FunctionTolerance', 1e-12));
    [m_conical_n, m_folded_n, v_conical_n, v_folded_n] = generateCandidatesGeneral(m_conical, m_folded, v_conical, v_folded, i, R, A, d, z_n(1), z_n(2), phi);    
end
function [m_conical_n_cand, m_folded_n_cand, v_conical_n_cand, v_folded_n_cand] = generateCandidatesGeneral(m_conical, m_folded, v_conical, v_folded, i, R, A, d, mz_n, vz_n, phi)
    m_folded_n_cand = (R^i)*[A + d*(2*i - 1); 0; mz_n];    
    v_folded_n_cand = (R^i)*[A + d*2*i; 0; vz_n];
    [m_conical_n_cand, ~] = throwTriangleConical_v2(m_folded(:, i), R*v_folded(:, i), m_folded_n_cand, m_conical(:, i), R*v_conical(:, i),   phi);
    [v_conical_n_cand, ~] = throwTriangleConical_v2(v_folded(:, i), m_folded_n_cand, v_folded_n_cand, v_conical(:, i), m_conical_n_cand,  phi);
end
function diff = helperGeneral(m_conical, m_folded, v_conical, v_folded, i, R, A, d, mz_n, vz_n, phi)    
    [m_conical_n_cand, m_folded_n_cand, v_conical_n_cand, v_folded_n_cand] = generateCandidatesGeneral(m_conical, m_folded, v_conical, v_folded, i, R, A, d, mz_n, vz_n, phi);
    diff = [norm(m_conical_n_cand - R*v_conical_n_cand) - norm(m_folded_n_cand - R*v_folded_n_cand),...
           norm(m_conical_n_cand - v_conical(:, i)) - norm(m_folded_n_cand - v_folded(:, i))];
end