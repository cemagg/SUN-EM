function [V] = FillVVector(rho_c_pls,rho_c_mns,EMag,theta_0,phi_0,dof2edge)
% FILLVVECTOR Fill the RHS vector.
%   At present, only a normally incident x-directed E field is implemented.
%   For others, eq 22 and 23 [RWG82] must still be implemented.

global NUM_DOFS EDGECONXELEMS ELL 

V = zeros(NUM_DOFS,1); % RHS matrix

E = zeros(NUM_DOFS,3); % Sampled incident electric field
% x-polarized field, unity magnitude
E(:,1) = EMag; % Special case
%E(:,2) = 1; % Special case

if theta_0 ~= 0 || phi_0 ~= 0 
    error('Only normally incident plane wave implemented at present ')
end

for mm = 1:NUM_DOFS
    pp_pls = EDGECONXELEMS(mm,1);
    pp_mns = EDGECONXELEMS(mm,2);
    V(mm) = dot(E(pp_pls,:)',rho_c_pls(mm,:))/2 + dot(E(pp_mns,:)',rho_c_mns(mm,:))/2;
    V(mm) = ELL(dof2edge(mm))*V(mm);
end


end

