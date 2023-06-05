% Returns the scalar cross product of two 2-element vectors 
% a x b = a(1)*b(2) - a(2)*b(1)
function c = cross2(a, b)
    c = a(1)*b(2) - a(2)*b(1);
end