% Returns 3x3 rotation matrix R that will rotate vectors by an angle theta
% (in rad) about an axis a (a column vector)
function R = getRotationMatrixAxisAngle(a, theta)
    a = [a(1); a(2); a(3)];
    a = a/norm(a);
    R = cos(theta)*eye(3) + (1 - cos(theta))*(a*a') + sin(theta)*vectorCross(a);
end