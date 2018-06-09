function [L_i,W_i] = GL_quad_rule(num_pts,unit_length)

% Gauss-Legendre quadrature rules for integration on a line, with integration limits (-1,1) as usual.
% Formulas are given in integer form for maximum accuracy, see Wikipedia entry. Checked against
% available real-valued data [p77, J Jin,  "Finite Elements in Electromagnetics", 2nd ed, Wiley 2002] to 9 digits.
% Presently implemented: 1 to 5 point rules, with degrees of precision
% 1,3,5,7 and 9 respectively (n-pt G-L quadrature has accuracy 2n-1).

% The optional flag unit_length, if present and true (1), returns the sampling points
% and weights scaled for integration limits (0,1), as is usual when working with
% simplex coordinates.

switch num_pts
    case 1 % p=1 accuracy
        L_i=0;
        W_i = 2;
    case 2 % p=3 accuracy
        L_i=[-1/sqrt(3)
            1/sqrt(3)];
        W_i = [1
            1];
    case 3 % p=5 accuracy
        L_i=[-sqrt(3/5)
            0
            sqrt(3/5)];
        W_i = [5/9
            8/9
            5/9];
    case 4 % p=7
        L_i=[-sqrt((3+2*sqrt(6/5))/7)
            -sqrt((3-2*sqrt(6/5))/7)
            sqrt((3-2*sqrt(6/5))/7)
            sqrt((3+2*sqrt(6/5))/7)];
        W_i = [(18-sqrt(30))/36
            (18+sqrt(30))/36
            (18+sqrt(30))/36
            (18-sqrt(30))/36];
    case 5 % p=9
        L_i=[-sqrt(5+2*sqrt(10/7))/3
            -sqrt(5-2*sqrt(10/7))/3
            0
            sqrt(5-2*sqrt(10/7))/3
            sqrt(5+2*sqrt(10/7))/3];
        W_i = [(322-13*sqrt(70))/900
            (322+13*sqrt(70))/900
            128/225
            (322+13*sqrt(70))/900
            (322-13*sqrt(70))/900];
    otherwise
        error('Unimplemented quadrature rule requested in tri_quad_rule')
end

if nargin==2 % scale to (0,1) integration limits
    if unit_length
        L_i = 0.5*L_i+0.5;
        W_i = 0.5*W_i;
    end
end
