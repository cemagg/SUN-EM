function PostProcMoM(I,EMag,dof2edge,eta_0,L,W,Xmesh,Ymesh,ProbType,quad_pts,sing)
%POSTPROCMOM Post process results

global NODE_COORD EDGES NUM_ELEMENTS NUM_EDGES NUM_DOFS ELL 

I_norm = eta_0*I/EMag; % Scale normalized current.

ell_vec = zeros(NUM_EDGES,3);
% post-process - 2D only
for iedge=1:NUM_EDGES
    ell_vec(iedge,:) = NODE_COORD(EDGES(iedge,2),:) - NODE_COORD(EDGES(iedge,1),:);
    ell_vec(iedge,:) = ell_vec(iedge,:)/ELL(iedge);
end

% Plot x directed currents on cut at x=L/2.
% Needs mesh with even number of elements in X direction; not checked!!
tol=max(ELL)/1e6;
ycounter = 0;
if rem(Xmesh,2)
    warning('Mesh increment in X direction must be even for plot AA')
else
    for idof=1:NUM_DOFS
        if abs(ell_vec(dof2edge(idof),2)-1) < tol && ...
                abs(NODE_COORD(EDGES(dof2edge(idof),1),1)-L/2) < tol
            ycounter = ycounter+1;
            %idof
            InormAA(ycounter) = I_norm(idof);
            y_vert(ycounter)=0.5*(NODE_COORD(EDGES(dof2edge(idof),1),2)+NODE_COORD(EDGES(dof2edge(idof),2),2));
        end
    end
end

xcounter = 0;
if rem(Ymesh,2)
    for idof=1:NUM_DOFS
        y_c = 0.5*(NODE_COORD(EDGES(dof2edge(idof),1),2)+NODE_COORD(EDGES(dof2edge(idof),2),2));
        if abs(ell_vec(dof2edge(idof),2)-1) < tol && ...
                abs(y_c-L/2) < tol
            xcounter = xcounter+1;
            %idof
            InormBB(xcounter) = I_norm(idof);
            x_hor(xcounter)=0.5*(NODE_COORD(EDGES(dof2edge(idof),1),1)+NODE_COORD(EDGES(dof2edge(idof),2),1));
        end
    end
else
    warning('Mesh increment in Y direction must be odd for plot BB')
end

figure
switch ProbType
    case 5
        ref_data_x_AA = L*linspace(1/10, 9/10, 5);
        ref_data_y_AA = [1.8 1.1 1.1 1.1 1.8]; % Still to fill in
        ref_data_x_BB = W*linspace(1/6, 5/6, 5);
        ref_data_y_BB = [0.8 1 1.1 1 0.8];
        max_y = 3;
    case 6
        ref_data_x_AA = L*linspace(1/14, 13/14, 7);
        ref_data_y_AA = [4.6 2.6 2.8 2.9 2.8 2.6 4.6];
        ref_data_x_BB = W*linspace(1/6, 5/6, 5);
        ref_data_y_BB = [1.7 2.5 2.9 2.5 1.7];
        max_y = 6;
end
if exist('InormAA','var') && exist('InormBB','var')
    plot(y_vert,abs(InormAA),'xb',ref_data_x_AA,ref_data_y_AA,'-ob',...
        x_hor,abs(InormBB),'+r',ref_data_x_BB,ref_data_y_BB,'--sr')
    if ~sing
    legend(['Cut AA: This code with ',num2str(NUM_ELEMENTS),' triangles, ',num2str(quad_pts),' quad pts '],...
        'Cut AA [Lit]',...
        'Cut BB this code, ditto',...
        'Cut BB [Lit]')
    else
        legend(['Cut AA: This code with ',num2str(NUM_ELEMENTS),' triangles, ',num2str(quad_pts),' quad pts  & sing.'],...
        'Cut AA [Lit]',...
        'Cut BB this code, ditto',...
        'Cut BB [Lit]')
    end
elseif exist('InormBB','var')
    plot(x_hor,abs(InormBB),'x',ref_data_x_BB,ref_data_y_BB,'o')
    legend(['Cut BB: this code with ',num2str(quad_pts),' quad pts'],'Cut BB [Lit]')
else
    plot(y_vert,abs(InormAA),'x',ref_data_x_AA,ref_data_y_AA,'o')
    legend(['Cut AA: this code with ',num2str(quad_pts),' quad pts'],'Cut AA [Lit]')
end
%title(['Dominant current component,',num2str(L),'\lambda square flat plate'])
axis([0 L 0 max_y])
xlabel('X/\lambda, Y/\lambda')
ylabel('|J_x/H_{inc}|')

plotfile = ['RWG_Fig_',num2str(ProbType),'_',num2str(NUM_DOFS),'dofs'];
print('-deps',plotfile)

end

