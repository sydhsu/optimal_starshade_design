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
%   nodes_planar - 3-D coordinates in unfolded/flat crease pattern
%   nodes_folded - 3-D coordinates in folded crease pattern
%   edges - array of node pairs that form an edge
%   triangulated - groups of 3 nodes after triangulating the entire pattern
%   tri_faces - groups of 3 nodes that form a triangular panel
%   quad_faces - groups of 4 nodes that form a quadrilateral panel
% ======================================================================= %

function [nodes_planar, nodes_folded, edges, triangulated, tri_faces, quad_faces] = flasher(N, n, h, A)
    beta    = 2*pi()/N;
    delta_z = 2*A*sin(beta/2)*tan(beta/2); % changes in height of mountain for a zero-thickness flasher
    d       = h/cos(beta/2);
    rot     = [ cos(beta), -sin(beta), 0;...
                sin(beta), cos(beta), 0;...
                0, 0, 1];

    m_nodes_planar = nan(3, n);
    m_nodes_folded = nan(3, n);
    v_nodes_planar = nan(3, n);
    v_nodes_folded = nan(3, n);

    m_nodes_planar(:, 1) = [A; 0; 0];
    m_nodes_folded(:, 1) = m_nodes_planar(:, 1);
    v_nodes_planar(:, 1) = m_nodes_planar(:, 1);
    v_nodes_folded(:, 1) = m_nodes_planar(:, 1);

    for i = 2:n
        v_nodes_folded(:, i) = (rot^(i - 1))*[A + 2*d*(i - 1); 0; 0];
    end

    [m_nodes_planar(:, 2), m_nodes_folded(:, 2), v_nodes_planar(:, 2)] = ...
        flasherSolveFirstPoint(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, delta_z);
    for i = 2:(n - 1)
        [m_nodes_planar(:, i + 1), m_nodes_folded(:, i + 1), v_nodes_planar(:, i + 1)] = ...
            flasherSolveGeneralPoint(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, i, rot, A, d, delta_z);
    end
    
    nodes_planar    = [];
    nodes_folded    = [];
    edges           = [];
    triangulated    = [];
    tri_faces       = [];
    quad_faces      = [];

    for j = 1:N
        nodes_planar    = [nodes_planar, (rot^(j - 1))*v_nodes_planar, (rot^(j - 1))*m_nodes_planar];
        nodes_folded    = [nodes_folded, (rot^(j - 1))*v_nodes_folded, (rot^(j - 1))*m_nodes_folded];
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
                m_id(2, j),     v_id(2, j)]; %;...      % minor mountain                    
%                 m_id(2, j),     v_id(2, j + 1)];    % diagonal
                            
        triangulated = [triangulated;...
                v_id(1, j),     m_id(2, j), v_id(2, j);...                    
                v_id(1, j),     m_id(2, j), v_id(1, j + 1);...
                v_id(1, j + 1), m_id(2, j), v_id(2, j + 1)]; % cuts quadrilateral panel in half

        tri_faces = [tri_faces;...
                v_id(1, j),     m_id(2, j), v_id(2, j);...                    
                v_id(1, j),     m_id(2, j), v_id(1, j + 1)];

        for i = 2:(n - 1)            
            edges = [edges;...
                    v_id(i, j),     v_id(i + 1, j);...  % major valley
                    m_id(i, j),     m_id(i + 1, j);...  % major mountain
                    m_id(i + 1, j), v_id(i, j + 1);...  % minor valley
                    m_id(i + 1, j), v_id(i + 1, j)];%;...  % minor mountain
%                     m_id(i, j),     v_id(i + 1, j);...  % diagonal                    
%                     m_id(i + 1, j), v_id(i + 1, j + 1)];% diagonal            
            triangulated = [triangulated;...
                    v_id(i, j),     m_id(i, j),     v_id(i + 1, j);...
                    m_id(i, j),     v_id(i + 1, j), m_id(i + 1, j);...
                    m_id(i, j),     m_id(i + 1, j), v_id(i, j + 1);...
                    v_id(i, j + 1), v_id(i + 1, j + 1), m_id(i + 1, j)];
            quad_faces = [quad_faces;...
                    v_id(i, j),     m_id(i, j),     m_id(i + 1, j),  v_id(i + 1, j);... 
                    m_id(i, j),     v_id(i-1,j+1),  v_id(i,j+1),     m_id(i+1, j)];
        end 
        edges = [edges;[m_id(n,j), v_id(n,j+1)]]; % adds outermost edge to close the outer polygon shape
    end
end

function [m_node_planar2, m_node_folded2, v_node_planar2] = flasherSolveFirstPoint(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, delta_z)
    fun = @(z2) helperFirst(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, z2);    
    z2  = fsolve(fun, delta_z, optimoptions('fsolve', 'Display', 'off'));    
    [m_node_planar2, m_node_folded2, v_node_planar2] = ...
        generateCandidatesFirst(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, z2);   
end

function [m_node_planar2, m_node_folded2, v_node_planar2] = generateCandidatesFirst(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, z2)
    m_node_folded2      = rot*[A + d; 0; z2];    
    [~, m_node_planar2] = throwTriangle(v_nodes_folded(:, 1), rot*v_nodes_folded(:, 1), m_node_folded2, v_nodes_planar(:, 1), rot*v_nodes_planar(:, 1));    
    [~, v_node_planar2] = throwTriangle(m_nodes_folded(:, 1), m_node_folded2, v_nodes_folded(:, 2), m_nodes_planar(:, 1), m_node_planar2);
end

function diff = helperFirst(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, z2)
    [m_node_planar2, m_node_folded2, v_node_planar2] = generateCandidatesFirst(m_nodes_planar, m_nodes_folded, v_nodes_planar, v_nodes_folded, rot, A, d, z2);    
    diff = norm(m_node_planar2 - rot*v_node_planar2) - norm(m_node_folded2 - rot*v_nodes_folded(:, 2));    
end

function [m_next, mp_next, v_next] = flasherSolveGeneralPoint(m, mp, v, vp, i, R, A, d, delta_z) 
    fun     = @(z_next) helperGeneral(m, mp, v, vp, i, R, A, d, z_next);
    options = optimoptions('fsolve', 'Display', 'off');
    z_prev  = mp(3, i);
    z_next  = fsolve(fun, z_prev + delta_z, options);
    [m_next, mp_next, v_next] = generateCandidatesGeneral(m, mp, v, vp, i, R, A, d, z_next);    
end

function [m_cand_next, mp_cand_next, v_cand_next] = generateCandidatesGeneral(m, mp, v, vp, i, R, A, d, z_next)
    mp_cand_next        = (R^i)*[A + d*(2*(i + 1) - 3); 0; z_next];
    [~, m_cand_next]    = throwTriangle(mp(:, i), R*vp(:, i), mp_cand_next, m(:, i), R*v(:, i));
    [~, v_cand_next]    = throwTriangle(vp(:, i), mp_cand_next, vp(:, i + 1), v(:, i), m_cand_next);
end

function diff = helperGeneral(m, mp, v, vp, i, R, A, d, z_next)    
    [m_cand_next, mp_cand_next, v_cand_next] = generateCandidatesGeneral(m, mp, v, vp, i, R, A, d, z_next);
    diff = norm(m_cand_next - R*v_cand_next) - norm(mp_cand_next - R*vp(:, i + 1));
end