% outputs the crease pattern (both unfolded and folded) for a
% thickness-accommodating flasher with N gores (or 2N major folds), n nodes
% per major fold line, h is the thickness, A is the radius of the hub, k is
% the fraction between 0 and 1 where the mountain node sits in the packaged
% configuration between the two valley nodes

function flasher_generate_cp_file(N, n, h, A, k, filename)
    fileID  = fopen(filename, 'w');

    beta    = 2*pi()/N;
    delta_z = 2*A*sin(beta/2)*tan(beta/2); % changes in height of mountain for a zero-thickness flasher
    R       = [ cos(beta), -sin(beta), 0;...
                sin(beta), cos(beta), 0;...
                0, 0, 1];

    m   = nan(3, n);
    mp  = nan(3, n);
    v   = nan(3, n);
    vp  = nan(3, n);

    m(:, 1)     = [A; 0; 0];
    mp(:, 1)    = m(:, 1);
    v(:, 1)     = m(:, 1);
    vp(:, 1)    = m(:, 1);

    for i = 2:n
        vp(:, i) = (R^(i - 1))*[A + (i - 1)*h; 0; 0];
    end

    [m(:, 2), mp(:, 2), v(:, 2)] = flasherSolveFirstPoint(m, mp, v, vp, R, A, k, h, delta_z);
    for i = 2:(n - 1)
        [m(:, i + 1), mp(:, i + 1), v(:, i + 1)] = flasherSolveGeneralPoint(m, mp, v, vp, i, R, A, k, h, delta_z);
    end
    
    nodes   = [];
    nodes_p = [];
    edges   = [];

    for j = 1:N
        nodes   = [nodes,   (R^(j - 1))*v,  (R^(j - 1))*m];
        nodes_p = [nodes_p, (R^(j - 1))*vp, (R^(j - 1))*mp];
    end
    
    i_grid = (ones(N, 1)*(1:n))';   % node id along the gore
    j_grid = ones(n, 1)*(1:N);      % node id across the circle
    
    v_id = n*(2*j_grid - 2) + i_grid;
    m_id = n*(2*j_grid - 1) + i_grid;  
    v_id = [v_id, v_id(:, 1)];
    m_id = [m_id, m_id(:, 1)];
    
    formatSpec = '%d %.1f %.1f %.1f %.1f\r\n';
    
    % the .cp format: (linetype: 1.Contour, 2.Mountain, 3.Valley) (start x) (start y) (end x) (end y)
           
    for j = 1:N        
        edges = [edges; 3, m_id(1, j), v_id(1, j + 1)];                
        for i = 1:(n - 1)
            edges = [edges;...
                    3, v_id(i, j),     v_id(i + 1, j);...  % major valley
                    2, m_id(i, j),     m_id(i + 1, j);...  % major mountain
                    3, m_id(i + 1, j), v_id(i, j + 1);...  % minor valley
                    2, m_id(i + 1, j), v_id(i + 1, j);...  % minor mountain
                    2, m_id(i, j),     v_id(i + 1, j);...  % diagonal mountain
                    3, m_id(i + 1, j), v_id(i + 1, j + 1)];% diagonal valley           
        end
    end
    
    for i = 1:size(edges, 1)                
        start_x = nodes(1, edges(i, 2));
        start_y = nodes(2, edges(i, 2));
        
        end_x = nodes(1, edges(i, 3));
        end_y = nodes(2, edges(i, 3));
        
        %display(sprintf(formatSpec, edges(i, 1), start_x, start_y, end_x, end_y));
        fprintf(fileID, formatSpec, edges(i, 1), start_x, start_y, end_x, end_y);
    end
    
    fclose(fileID);
end