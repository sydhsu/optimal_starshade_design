function [edgeNodes, edgeEdges] = opticalShieldEdgework(nodes_planar, N, n, t_p_percent, R_truss_deployed)
    % step 2
    i_v_endm1 = n - 1 + 2*n - 1;
    i_v_end   = n + 2*n - 1;

    i_m_endm1 = 2*n - 1;
    i_m_end   = 2*n;

    v_endm1    = nodes_planar(1:2, i_v_endm1);
    v_end      = nodes_planar(1:2, i_v_end);
    m_endm1    = nodes_planar(1:2, i_m_endm1);
    m_end      = nodes_planar(1:2, i_m_end);
    
    t_p = t_p_percent*norm(m_end - m_endm1);

    % step 3
    p   = t_p*(m_end - m_endm1)/norm(m_end - m_endm1) + m_endm1;

    rot = [cos(pi()/N), -sin(pi()/N); sin(pi()/N) cos(pi()/N)];
    A = [(v_end - v_endm1)/norm(v_end - v_endm1), -rot*p/norm(p)];
    x = A\(-v_endm1);
    q = x(1)*(v_end - v_endm1)/norm(v_end - v_endm1) + v_endm1;

    dt = pi()/N;
    t = (0:dt:(2*pi() - dt)) + atan2(p(2), p(1));
    X_truss_deployed = R_truss_deployed*cos(t);
    Y_truss_deployed = R_truss_deployed*sin(t);

    truss_m = [X_truss_deployed(1); Y_truss_deployed(1)];
    truss_v = [X_truss_deployed(2); Y_truss_deployed(2)];

    % step 4
    x1 = dot((p - v_endm1), (v_end - v_endm1))/norm(v_end - v_endm1);
    s = x1*(v_end - v_endm1)/norm(v_end - v_endm1) + v_endm1;

    angle_1 = angleBetweenVectors(q - p, s - p);

    alpha = angleBetweenVectors(v_endm1 - v_end, p - q);

    % step 5
    t = (pi() - alpha)/2;
    v_hat = [cos(t), sin(t); -sin(t), cos(t)]*truss_v/norm(truss_v);
    w_hat = [cos(angle_1), sin(angle_1); -sin(angle_1), cos(angle_1)]*(q - p)/norm(q - p);
    A = [v_hat, -w_hat];
    x = A\(-q + p);
    a1 = x(1)*v_hat + q;

    beta = angleBetweenVectors(q - a1, p - a1);

    q_hat = q/norm(q);
    a1_r = 2*dot(a1 - q, q_hat)*q_hat - (a1 - q) + q;

    % step 6 and 7
    A = [(truss_v - truss_m)/norm(truss_v - truss_m), -(a1 - p)/norm(a1 - p)];
    x = A\(-truss_m + p);

    O2 = x(2)*(a1 - p)/norm(a1 - p) + p;

    angle_2 = (1/2)*angleBetweenVectors(truss_v -  truss_m, a1 - p);

    v_hat = [cos(beta), -sin(beta); sin(beta), cos(beta)]*(a1 - a1_r)/norm(a1 - a1_r);
    w_hat = [cos(angle_2), sin(angle_2); -sin(angle_2), cos(angle_2)]*(a1 - O2)/norm(a1 - O2);
    A = [v_hat -w_hat];
    x = A\(-a1 + O2);
    a2 = x(1)*v_hat + a1;

    a2_r = 2*dot(a2 - q, q_hat)*q_hat - (a2 - q) + q;

    gamma = angleBetweenVectors(a2_r - a2, a1 - a2);

    % step 8
    v_hat = ([cos(pi() - gamma), -sin(pi() - gamma); sin(pi() - gamma), cos(pi() - gamma)]*(O2 - a2))/norm(O2 - a2);
    A = [v_hat, -(truss_v - truss_m)/norm(truss_v - truss_m)];
    x = A\(-a2 + truss_m);
    a3 = x(1)*v_hat + a2;

    a3_r = 2*dot(a3 - q, q_hat)*q_hat - (a3 - q) + q;

    % step 9 and 10
    angle_3 = (1/2)*angleBetweenVectors(truss_v - truss_m, -p);
    v_hat = [cos(angle_3), -sin(angle_3); sin(angle_3), cos(angle_3)]*(truss_v - truss_m)/norm(truss_v - truss_m);
    A = [v_hat, -(O2 - a2)/norm(O2 - a2)];
    x = A\(-truss_m + a2);
    b1 = x(1)*v_hat + truss_m;

    p_hat = p/norm(p);
    b1_r = 2*dot(b1 - p, p_hat)*p_hat - (b1 - p) + p;
    
    silentFlag = true;
    
    if not(silentFlag)
        disp('Step 3 check')
        disp(['Angle between p and q: ' num2str(angleBetweenVectors(p, q))])
        disp(['Angle between t_m and t_v: ' num2str(angleBetweenVectors(truss_m, truss_v))])
        disp(['pi/N: ' num2str(pi()/N)])

        disp('Step 4 check')
        disp(['Angle between p, s and v_n-1, v_n: ' num2str(angleBetweenVectors(p - s, v_end - v_endm1))])

        disp('Step 5 check')
        disp(['Angle a_1, p, q: ' num2str(angleBetweenVectors(a1 - p, q - p))])
        disp(['Angle q, p, s: ' num2str(angleBetweenVectors(q - p, s - p))])
        disp(['Angle a1, q, a1r: ' num2str(angleBetweenVectors(a1 - q, a1_r - q)) '; should be: ' num2str(pi() - alpha)])
        disp(['Angle a1, q, truss_v: ' num2str(angleBetweenVectors(a1 - q, truss_v - q)) '; should be: ' num2str((pi() - alpha)/2)])

        disp('Step 6 and 7 check')
        disp(['1/2 of angle between tv - tm, a1 - p: ' num2str(0.5*angleBetweenVectors(truss_v - truss_m, a1 - p))])
        disp(['Angle between a2 - O2, a1 - p: ' num2str(angleBetweenVectors(a2 - O2, a1 - p))])
        disp(['Angle a1r, a1, a2: ' num2str(angleBetweenVectors(a2 - a1, a1_r - a1)) '; should be: ' num2str(pi() - beta)])

        disp('Step 8 check')
        disp(['Angle b1, a2, a3: ' num2str(angleBetweenVectors(b1 - a2, a3 - a2)) '; should be: ' num2str(pi() - gamma)])

        disp('Step 9 check')
        disp(['Angle b1, tm, tv: ' num2str(angleBetweenVectors(truss_v - truss_m, b1 - truss_m))])
        disp(['Angle b1, tm, p: ' num2str(angleBetweenVectors(b1 - truss_m, p - truss_m))])
    end
               %1         2  3  4   5   6    7      8     9     10       11     12 13  14     15         
    edgeNodes = [v_endm1, s, q, a1, a2, a3, a1_r, a2_r, a3_r, truss_v, m_endm1, p, b1, b1_r, truss_m];
    edgeEdges = [1, 2; 2, 3; 3, 4; 3, 7; 4, 7; 4, 5; 7, 8; 5, 8; 5, 6; 8, 9; 6, 9; 6, 10; 9, 10; 11, 12; 12, 13; 12, 14; 13, 14; 13, 15; 14, 15; 6, 15; 5, 13; 4, 12; 3, 12; 2, 12];

end