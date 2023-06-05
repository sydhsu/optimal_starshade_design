% returns a 3x3 matrix ax that has the following form for a 1x3 vector a

function ax = vectorCross(a)
    ax = [  0, -a(3), a(2);...
            a(3), 0, -a(1);...
            -a(2), a(1), 0];
end