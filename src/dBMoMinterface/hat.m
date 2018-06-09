function [hat_a]=hat(a)
% HAT returns hat(a), the unit vector in the direction of vector a.
% If a is 0 to working precision, then the zero vector of the same length as a is returned.
if abs(a) < eps
    hat_a = 0*a;
else
    hat_a = a/norm(a);
end
    