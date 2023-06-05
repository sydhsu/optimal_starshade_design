% returns the 2d angle (in radians) between the lines [P0,P1] and [P0,P2]

function angle = angleBetweenPoints(P0,P1,P2)
    angle = angleBetweenVectors(P1-P0, P2-P0);
end