% Takes the triangle defined by p1, p2, p3 and positions this triangle such
% that in its new position p1p, p2p, p3p, p3p has z-value z0

% Returns two possible positions of p3p
% Returned values are nans if there are no solutions

function [p3p_plus, p3p_minus] = positionTriangleByZ(p1, p2, p3, p1p, p2p, z0)
    axis    = (p2 - p1)/norm(p2 - p1);
    q       = dot(p3 - p1, p2 - p1)/norm(p2 - p1)*axis + p1;
    r       = norm(p3 - q);
    
    axisp   = (p2p - p1p)/norm(p2p - p1p);
    qp      = dot(p3 - p1, p2 - p1)/norm(p2 - p1)*axisp + p1p;
    
    axp     = axisp(1);
    ayp     = axisp(2);
    azp     = axisp(3);
    
    vzp     = z0 - qp(3);
    
    vxp     = roots([(1 + (axp/ayp)^2), 2*vzp*axp*azp/(ayp^2), ((vzp*azp/ayp)^2 + vzp^2 - r^2)]);
    
    if not(isreal(vxp))
        disp('No solutions')
        vxp = nan(1, 2);    
    end
    
    vyp     = -(vxp*axp + vzp*azp)/ayp;
    
    p3p_plus    = [vxp(1); vyp(1); vzp] + qp;
    p3p_minus   = [vxp(2); vyp(2); vzp] + qp;
    
%     disp(['Plus error to p1: ' num2str( (norm(p3p_plus - p1p) - norm(p3 - p1))/norm(p3 - p1)*100 ) ])
%     disp(['Plus error to p2: ' num2str( (norm(p3p_plus - p2p) - norm(p3 - p2))/norm(p3 - p2)*100 ) ])
%     
%     disp(['Minus error to p1: ' num2str( (norm(p3p_minus - p1p) - norm(p3 - p1))/norm(p3 - p1)*100 ) ])
%     disp(['Minus error to p2: ' num2str( (norm(p3p_minus - p2p) - norm(p3 - p2))/norm(p3 - p2)*100 ) ])
end