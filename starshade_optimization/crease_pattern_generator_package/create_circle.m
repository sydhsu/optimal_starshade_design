function [x,y] = create_circle(R,x0,y0,n)    
    dt = 2*pi()/n;
    t = 0:dt:(2*pi()-dt);
    
    x = R*cos(t) + x0;
    y = R*sin(t) + y0;
end