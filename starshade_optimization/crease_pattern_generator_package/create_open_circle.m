function [x,y] = create_open_circle(R,x0,y0,theta1,theta2,n)        
    t = linspace(theta1, theta2, n);
    
    x = R*cos(t) + x0;
    y = R*sin(t) + y0;
end