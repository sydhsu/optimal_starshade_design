close all; clear all; clc;

for i = 1:100
    %p1 = zeros(3, 1);
    %p2 = p1 + [rand(1); 0; 0];
    %p3 = p1 + [rand(1); rand(1); 0];

    p1 = rand(3, 1);
    p2 = rand(3, 1);
    p3 = rand(3, 1);

    r1 = 0.8*rand(1);
    r2 = 0.8*rand(1);
    r3 = 1.1*rand(1);

    [p_plus, p_minus] = intersectionOfThreeSpheres(p1, p2, p3, r1, r2, r3);

    disp(['Plus to 1: ' num2str(norm(p_plus - p1) - r1)])
    disp(['Plus to 2: ' num2str(norm(p_plus - p2) - r2)])
    disp(['Plus to 3: ' num2str(norm(p_plus - p3) - r3)])

    disp(['Minus to 1: ' num2str(norm(p_minus - p1) - r1)])
    disp(['Minus to 2: ' num2str(norm(p_minus - p2) - r2)])
    disp(['Minus to 3: ' num2str(norm(p_minus - p3) - r3)])
end
