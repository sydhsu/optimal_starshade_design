% Computes the rigid body translation and rotation between 3 points
% Non-rigid deformations will cause this to fail spectacularly

% P contains 3 position column vectors
% P_prime contains 3 position column vectors

% Returns rotation matrix C and translation L such that P = C*P_prime + L
% the error is an RMS distance error = mean(norm(P - (C*P_prime + L))

function [C, L, err] = matchPlanes(P,P_prime)
    
    % find the relative in-plane vectors between the position vectors
    V = P(:,[2 3]) - P(:,[1 1]);
    V_prime = P_prime(:,[2 3]) - P_prime(:,[1 1]);
    
    % find the surface normals
    n = cross(V(:,1),V(:,2));
    n = n/norm(n);
    n_prime = cross(V_prime(:,1),V_prime(:,2));
    n_prime = n_prime/norm(n_prime);
    
    % compute the rotation matrix C1 so that V and C1*V_prime lie in the same
    % plane
    a = cross(n,n_prime);
    a = a/norm(a);
    cos_phi = dot(n,n_prime)/norm(n)/norm(n_prime);
    sin_phi = det([n'; n_prime'; a'])/norm(n)/norm(n_prime);
    C1 = cos_phi*eye(3) + (1 - cos_phi)*(a*a') - sin_phi*[0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];
              
    % now rotate about the common surface normal so that V and V_prime
    % align    
    a = n;
    cos_phi = dot(V(:,1),C1*V_prime(:,1))/norm(V(:,1))/norm(C1*V_prime(:,1));
    sin_phi = det([V(:,1)'; (C1*V_prime(:,1))'; a'])/norm(V(:,1))/norm(C1*V_prime(:,1));
    C2 = cos_phi*eye(3) + (1 - cos_phi)*(a*a') - sin_phi*[0, -a(3), a(2); a(3), 0, -a(1); -a(2), a(1), 0];    
           
    C = C2*C1;
    
    L = mean(P - C*P_prime,2);    
    P_est = C*P_prime + L*ones(1,3);
    
    err_2 = 0;
    for i = 1:3
        err_2 = err_2 + norm(P(:,i) - P_est(:,i))^2;
    end
    err = sqrt(err_2/3);
    
end