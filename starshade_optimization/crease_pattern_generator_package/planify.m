% takes a set of 3D points that defines a chain of triangular faces in 3D
% and outputs a set of points on a plane that defines an isometric chain of
% triangular faces on a plane. The planar chain of triangles, when folded
% along the edges of the triangles, can reproduce the 3D configuration of
% the triangles

function points2d = planify(points3d)
    n = size(points3d, 2);
    
    points2d = nan(3, n);
    
    points2d(:, 1) = [0; 0; 0];
    points2d(:, 2) = [norm(points3d(:, 2) - points3d(:, 1)); 0; 0];
    
    for i = 1:(n - 2)        
        b1 = points3d(:, i + 1) - points3d(:, i);
        b1 = b1/norm(b1);        
        b2 = (points3d(:, i + 2) - points3d(:, i)) - dot((points3d(:, i + 2) - points3d(:, i)), b1)*b1;
        b2 = b2/norm(b2);        
        
        a1 = (points2d(:, i + 1) - points2d(:, i));
        a1 = a1/norm(a1);
        a2 = [0 -1 0; 1 0 0; 0 0 1]*a1;        
        
        s = 1;
        
        if i > 1
            s = -sign(dot(points2d(:, i - 1) - points2d(:, i), a2));
        end
        
        points2d(:, i + 2) =    dot(points3d(:, i + 2) - points3d(:, i), b1)*a1 + ...
                                s*dot(points3d(:, i + 2) - points3d(:, i), b2)*a2 + ...
                                points2d(:, i);
    end
end