function [lambda]=simplex_area(q,a,b,c)
% SIMPLEX_AREA returns the simplex (area) coordinates of the point q
% with coordinates (q1,q2,q3) for a triangle with vertices a,b and
% c. If q is not in the plane of the triangle, the routine returns the
% simplex coordinates of the projected point in the plane.

% The algorithm (but not code) comes from
% [Press] Press et al, "Numerical Recipes: the Art of Scientific
% Computing", 3rd ed, CUP 2007, p.1116.

% alpha is the simplex coordinate defined from the edge opposite node a; 
% beta is the simplex coordinate defined from the edge opposite node b; 
% gamma (not computed explicitly) is the last, defined from the edge
% opposite node c.
% These correspond to the widespread usage lambda_1,lambda_2,lambda_3 for nodes 1,2,
% and 3 in computational EM, which is what the code returns.

ap=a-c;
bp=b-c;
qp=q-c;

denom = dot(ap,ap)*dot(bp,bp)-dot(ap,bp)^2;
alpha = (dot(bp,bp)*dot(ap,qp)-dot(ap,bp)*dot(bp,qp))/denom;
beta  = (dot(ap,ap)*dot(bp,qp)-dot(ap,bp)*dot(ap,qp))/denom;

lambda=[alpha; beta; 1-alpha-beta];
