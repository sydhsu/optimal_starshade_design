close all; clear all; clc;

for i = 1:100

    O = rand(3, 1);
    a1 = rand(3, 1);
    a2 = rand(3, 1);
    L = rand(1);

    theta1 = pi()/2*rand(1);
    theta2 = 0.1*rand(1) + (acos(dot(a1, a2)/norm(a1)/norm(a2)) - theta1);

    [p_plus, p_minus] = intersectionOfTwoConesWithTheSameApex(O, a1, theta1, a2, theta2, L);

    disp(['Plus theta 1: '    num2str(theta1 - acos(dot((p_plus -  O), a1)/norm(p_plus -  O)/norm(a1)))])
    disp(['Plus theta 2: '    num2str(theta2 - acos(dot((p_plus -  O), a2)/norm(p_plus -  O)/norm(a2)))])
    disp(['Plus L: ' num2str(norm(p_plus - O) - L)])
    disp(['Minus theta 1: '   num2str(theta1 - acos(dot((p_minus - O), a1)/norm(p_minus - O)/norm(a1)))])
    disp(['Minus theta 2: '   num2str(theta2 - acos(dot((p_minus - O), a2)/norm(p_minus - O)/norm(a2)))])
    disp(['Minus L: ' num2str(norm(p_plus - O) - L)])

end