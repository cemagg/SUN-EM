% MOM3D_DEMO uses RWG basis functions to solve for the surface current
% induced on a PEC scatterer. The code uses the orignal quasi-Galerkin approach
% [RWG82] below. See [RWG82] SM Rao, DR Wilton & AW Glisson
% "Electromagnetic scattering by surfaces of arbitrary shape", IEEE
% Transactions on Antennas Propagation, May 1982. Vol AP-3, No 3, pp.
% 409-418.
%
% The present version of the code handles only a regular, triangular grid,
% generated from a rectangular block mesh, with an x-directed uniform plane
% wave normally incident. Two standard tests in [RWG82] can be run, and
% results are graphically compared to that reference.
%
% Mesh density is controlled by the 3rd and 4th arguments passed to
% function trimesh, very similar to that used for the 2D FEM solver. (The
% one called here also returns z coordinates, all set to z=0). The mesher
% is specialized, but works with both odd and even increments.
%
% The impedance matrix is using by edge (dof) and by face (element). 
% (This is to illustrate the different approaches; in practice, one or the other
% should be used). The latter approach re-uses many integrals, although the additional
% coding overhead of this approach results in less-speed up than indicated
% in the literature (at least, in this implementation).
%
% The singularity is handled either using quadrature, in an approximate fashion, 
% or using a modern singularity cancellation scheme.

% Author: DB Davidson, Stellenbosch Univ, June 2010. Inputs from D Ludick
% gratefully acknowledged.

clear all;
close all;

global ELEMENTS NODE_COORD EDGES ELEMENT_EDGES NUM_NODES NUM_ELEMENTS NUM_EDGES LOCALEDGENODES ...
    NUM_DOFS EDGECONXELEMS DOFLOCALNUM ELL LOCALVERTEX ELEMENT_PLS_MNS 

% Set up basic constants
eps_0 = 8.854e-12;
mu_0=4*pi*1e-7;
eta_0 = sqrt(mu_0/eps_0);

c = 1/sqrt(eps_0*mu_0); % Speed of light
freq = c; % Frequency in Hz
omega = 2*pi*freq;
lambda = c/freq; % Wavelength in m
k  = 2*pi/lambda; % Wavenumber in rad/m
EMag = 10; % Magnitude of incident field in V/m (not set to one to test code)
theta_0 = 0; % Angle of arrival of plane wav6e [degrees]
phi_0 = 0;

% Caution - do not change following three lines, hard-wired into other
% parts of code. Uses the convention of [Fig 4, RWG82]. Note that this is NOT
% the same as the standard FEM convention.
LOCALEDGENODES(1,:) = [2 3];
LOCALEDGENODES(2,:) = [1 3];
LOCALEDGENODES(3,:) = [1 2];
% LOCALVERTEX(j) Following gives the local vertex associated with edge j
% See [Fig 4, RWG82].
LOCALVERTEX(1) = 1;
LOCALVERTEX(2) = 2;
LOCALVERTEX(3) = 3;

% Set up geometry - initially, square plate lying in z=0 plane

% Mesh plate
ProbType = input('Enter problem to run: 5 (Fig 5, [RWG82}); 6 (Fig 6 - default)');
if isempty(ProbType)
    ProbType = 6;
end

switch ProbType
    case 5
        L=0.15*lambda; % [Fig 5, RWG82]
        W=L;
        Xmesh=6; % 6 for good results
        Ymesh=5; % 5 ditto
    case 6
        L=1*lambda;    % [Fig 6, RWG82]
        W=L;
        Xmesh=6; % 6
        Ymesh=7; % 7
    otherwise
        error('Unknown problem.')
end

% % Mesh plate
% Mult = input('Enter mesh blocks multiplier [default 1]: ');
% if isempty(Mult)
%     Mult = 1;
% end
% Xmesh=6*Mult; % 6
% Ymesh=7*Mult; % 7

sing = input('Use singularity integral scheme (T/F)? [default F]: ');
if isempty(sing)
    sing = 0;
end
quad_pts = 6;

tstart=tic;
[x_nodes,y_nodes,z_nodes] = trimesh3D(L,W,Xmesh,Ymesh); % generate triangular mesh
triplot(ELEMENTS,NODE_COORD(:,1),NODE_COORD(:,2));
axis([0 L 0 W]);
axis square
for inode = 1:NUM_NODES
    text(NODE_COORD(inode,1),NODE_COORD(inode,2),num2str(inode))
end
r_c = zeros(NUM_ELEMENTS,3); % Centres of elements in global system
for ielem = 1:NUM_ELEMENTS
    r_c(ielem,1) = mean(NODE_COORD(ELEMENTS(ielem,1:3),1));
    r_c(ielem,2) = mean(NODE_COORD(ELEMENTS(ielem,1:3),2));
    r_c(ielem,3) = mean(NODE_COORD(ELEMENTS(ielem,1:3),3)); % General
    text(r_c(ielem,1),r_c(ielem,2),num2str(ielem)); % Assumes plate is in z=0 plane.
end
print -deps tri_mesh

edgemake_MoM; % General
dof_int = outside_edge(L,W); % Assumes plate is in z=0 plane
[dof_RWG,dof2edge] = renumber_RWG(dof_int); % General
edge_conx_elem(dof_RWG); % General
find_local_dofs(dof_RWG); % General
figure
triplot(ELEMENTS,NODE_COORD(:,1),NODE_COORD(:,2));
axis([0 L 0 W]);
axis square
for idof = 1:NUM_DOFS
    edge_c = (NODE_COORD(EDGES(dof2edge(idof),1),:)+NODE_COORD(EDGES(dof2edge(idof),2),:))/2;
    text(edge_c(1),edge_c(2),num2str(idof))
end

% Compute vectors to centre points of triangles.
[rho_c_pls,rho_c_mns,x_grid_pls,y_grid_pls,x_grid_mns,y_grid_mns] = ComputeRho_c(r_c);

% Plot the above vectors (for debugging only).
hold % This plotting stub assumes z=0
quiver(x_grid_pls,y_grid_pls,rho_c_pls(:,1),rho_c_pls(:,2))
axis([0 L 0 W]);
quiver(x_grid_mns',y_grid_mns',rho_c_mns(:,1),rho_c_mns(:,2),'r')
TimePreProcessing = toc(tstart)

% Fill the impedance matrix
tstart1=tic;
[Z] = FillZMatrixByEdge(omega,eps_0,mu_0,k,r_c,rho_c_pls,rho_c_mns,quad_pts,sing,dof2edge);
TimeMatrixFillByEdge = toc(tstart1)

tstart1a=tic;
[Z_1] = FillZMatrixByFace(omega,eps_0,mu_0,k,r_c,rho_c_pls,rho_c_mns,quad_pts,sing,dof2edge,dof_RWG);
% report maximum difference
MatrixDiff=max(max(Z-Z_1))
TimeMatrixFillByFace = toc(tstart1a)


% Fill the RHS vector
[V] = FillVVector(rho_c_pls,rho_c_mns,EMag,theta_0,phi_0,dof2edge);

tstart2=tic;
I = Z\V;
TimeMatrixSolve = toc(tstart2)
I_1 =Z_1\V;
CurrDiff=max(max(I-I_1))

PostProcMoM(I,EMag,dof2edge,eta_0,L,W,Xmesh,Ymesh,ProbType,quad_pts,sing);
TimeOverall = toc(tstart)

savefile = ['test_',num2str(NUM_DOFS),'sing_',num2str(sing)];
save(savefile)

