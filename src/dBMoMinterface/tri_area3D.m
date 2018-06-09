function [area]=tri_area3D(a,b,c)
% TRI_AREA_3d returns the unsigned area of the triangle located in 3D space
% with vertices a, b, c. 

% The algorithm (but not code) comes from
% [Press] Press et al, "Numerical Recipes: the Art of Scientific
% Computing", 3rd ed, CUP 2007, p.1115.
vec_area = 0.5 *cross((a-c),(b-c));
area = norm(vec_area);