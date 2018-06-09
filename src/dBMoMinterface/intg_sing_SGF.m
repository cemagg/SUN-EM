function [Int,Int_xi1,Int_xi2,Int_xi3] = intg_sing_SGF(k,r,r1,r2,r3,num_pts_rad,num_pts_trans)
% This function evaulates the integrals I^pq, I_xi^pq or I_eta^pq at
% wavenumber k (rad/m) observation point r (all distances in m) over the triangle with nodes r1, r2 and r3 (and
% r(1),r(2),r(3) etc, corresponding to x,y and z coordinates). Point r may or
% may not lie in the plane of the triangle.
% The number of samples is determined by the last two arguments, which
% determine the radial (j) and transverse (i) samples respectively .
% These integrals involve the dynamic scalar Green's function in a region containing the singularity,
% as in Rao, Wilton and Glisson, "Electromagnetic scattering by surfaces of
% arbitrary shape", IEEE T-AP, Vol AP-30 No 3, May 1982, pp.409-418.
% The formulation uses a transform to map the triangular region into
% a quadrilateral, with the singularity softened and integrable using
% standard quadrature routines. Nested Gauss-Legendre quadrature is used,
% with num_pts_i the number of points in the u (mapped x) direction, and
% num_pts_j the number of points in the y direction.

% At present, the arcsinh rule has been implemented as the "singularity softener": see
% [K&W] Khayat & Wilton, "Numerical Evaluation of Singular and Near-Singular
% Potential Integrals", IEEE T-AP, Vol 53, Nov 10, Oct 2005, pp. 3180-3190.
% Other schemes, such as the original Duffy scheme, are discussed in
% Khayat, Wilton & Fink, "An Improved Transformation and Optimized Sampling Scheme
% for the Numerical Evaluation of Singular and Near-Singular Potentials",
% IEEE AWPL Vol 7 2008, pp. 377-380.

% Some of the computational geometry algorithms (but not code) come from
% [Press] Press et al, "Numerical Recipes: the Art of Scientific Computing", 3rd ed, CUP 2007.

% Author: DB Davidson, Stellenbosch University, April 2009.
%         Revised 22 Dec 2009 to return four integrals.

% WARNING: although the theory is valid for projected points falling
% outside the triangle, a quick test has showed that this may not be so with
% this implementation. This has not been investigated further.

num_pts_j=num_pts_rad;
num_pts_i=num_pts_trans;

% Zero integral
I(1:3) = 0;
I1(1:3) = 0;
I2(1:3) = 0;
I3(1:3) = 0;

% For initial testing, set basis function constant
%Lambda(1:num_pts_j,1:num_pts_i) = 1; % Now hard-coded below

% First, find projection of observation point r0 into plane of triangle:

N = cross((r2-r1),(r3-r1)); % [Press, eq.(21.3.17),p.1115], with a=r1,b=r2 and c=r3.
r0 = r+dot((r1-r),N)*N/(norm(N))^2; % [Ibid, eq.(21.3.18),p.1115], with p=r and p'=r0

% Find height of point above (below) plane of triangle. Construct a vector
% from projection of observation point to observation point and then find
% the projection vector onto the vector normal to the triangle.
r_vec = r-r0;
z=dot(hat(N),r_vec);

[xi0]=simplex_area(r0,r1,r2,r3); % simplex coordinates of projected observation point in original triangle.

xi_pij = zeros(3,1);
for tt=1:3 % Loop over sub-triangles 1,2 and 3
    % Compute geometry for relevant sub-triangle [K&W eq.4]
    r1_p = r0;
    switch tt % sub-triangle number
        case 1
            r2_p =r2;
            r3_p =r3;
            T = [xi0(1) 0 0; xi0(2) 1 0; xi0(3) 0 1]; %[K&W eq.17]
        case 2
            r2_p =r3;
            r3_p =r1;
            T = [xi0(1) 0 1; xi0(2) 0 0; xi0(3) 1 0]; %[K&W eq.17]
        case 3
            r2_p =r1;
            r3_p =r2;
            T = [xi0(1) 1 0; xi0(2) 0 1; xi0(3) 0 0]; %[K&W eq.17]
    end
    ell1_p=r3_p - r2_p;
    ell2_p=r1_p - r3_p;
    ell3_p=r2_p - r1_p;
    hat_n_p = hat(cross(ell1_p,ell2_p));
    A_p = dot(hat_n_p,cross(ell1_p,ell2_p))/2; % Expression for A_p corrected
    h1_p = 2*A_p/norm(ell1_p)^2*cross(ell1_p,hat_n_p);

    % Approximation of integral using nested Gauss-Legendre quadrature

    [xi_i,w_i]=GL_quad_rule(num_pts_i,1);
    [xi_j,w_j]=GL_quad_rule(num_pts_j,1);
    for jj= 1:num_pts_j % Note that [K&W eq.10] have the summation operators in the wrong sequence - j is the outer loop
        xi_pij(1)=xi_j(jj); % [K&W eq.13]
        y_pj =norm(h1_p)*(1-xi_j(jj)); %  [K&W eq.13]; h1_p here is the magnitude of the vector
        x_Lj =  dot(hat_n_p,cross(hat(h1_p),ell2_p))*(1-xi_pij(1)); %[K&W eq.6];
        x_Uj = -dot(hat_n_p,cross(hat(h1_p),ell3_p))*(1-xi_pij(1)); %[K&W eq.7];
        u_Lj = asinh(x_Lj/sqrt(y_pj^2+z^2)); %[K&W eq.9];
        u_Uj = asinh(x_Uj/sqrt(y_pj^2+z^2));
        for ii =1:num_pts_i
            u_ij = u_Lj*(1-xi_i(ii))+u_Uj*xi_i(ii); %[K&W eq.14];
            x_pij = sqrt(y_pj^2+z^2) *sinh(u_ij); %[K&W eq.15];
            xi_pij(3)= dot(hat_n_p,cross(ell3_p,(hat(h1_p)*y_pj- hat(ell1_p)*x_pij)))/(2*A_p); %[K&W eq.16];
            xi_pij(2)= 1 - xi_pij(3) - xi_pij(1);
            xi_ij = T*xi_pij; % [K&W, eq.17]; Original triangle area coordinates
            R_ij = norm(r-r1*xi_ij(1)-r2*xi_ij(2)-r3*xi_ij(3)); % [K&W, eq.18]
            % Finally, compute the integral in [K&W, eq.10]. Note that this
            % integral is computed first over sub-triangle 1, then over
            % sub-triangle 2 etc; these are recorded separately to aid
            % debugging, but could just be accumulated as one running total.
            %I(tt) = I(tt) + w_i(ii)*w_j(jj)*norm(h1_p)*(u_Uj-u_Lj)*Lambda(jj,ii)*exp(-j*k*R_ij);
            I(tt)  = I(tt)  + w_i(ii)*w_j(jj)*norm(h1_p)*(u_Uj-u_Lj)*exp(-j*k*R_ij);
            % Error may lie in next lines - this intergral is evaluated in
            % Cartesian space, but xi_ij are simplex coordinates????
            I1(tt) = I1(tt) + w_i(ii)*w_j(jj)*norm(h1_p)*(u_Uj-u_Lj)*xi_ij(1)*exp(-j*k*R_ij);
            I2(tt) = I2(tt) + w_i(ii)*w_j(jj)*norm(h1_p)*(u_Uj-u_Lj)*xi_ij(2)*exp(-j*k*R_ij);
            %Next line is redundant - testing only
            %I3(tt) = I3(tt) + w_i(ii)*w_j(jj)*norm(h1_p)*(u_Uj-u_Lj)*xi_ij(3)*exp(-j*k*R_ij);
        end
    end
end
% Add sub-triangle results and find last integral
Int     =sum(I);
Int_xi1 =sum(I1);
Int_xi2 =sum(I2);
Int_xi3 =Int-Int_xi1-Int_xi2;
%Int_xi3 =sum(I3);

