function [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = Int_pq(p,q,r_cp,k,quad_pts,sing)
% INT_PQ returns the four integrals for faces (=elements) p and q,
% corresponding to [eqs.34a-34d,RWG82]. Note that these integrals are
% evaluated in normalized coordinates, and are not scaled by the area.
% The overall matrix element is assembled from these elsewhere.
% Note also that p is the field point, q the source point. The integrals are
% performed over triangle q.
% Note further that xi, eta and zeta are equivalent to
% lambda_1, lambda_2, and lambda_3 respectively.
% If flag sing is set, then the singular terms are evaluated using a
% special integration rule.

% Author: D B Davidson, Dec 2009.
% Corrections for singular intergral evaluation: 1 June 2010 DBD.

global ELEMENTS NODE_COORD

qnodes = ELEMENTS(q,:);
n1 = NODE_COORD(qnodes(1),:);
n2 = NODE_COORD(qnodes(2),:);
n3 = NODE_COORD(qnodes(3),:);
area = tri_area3D(n1,n2,n3);

Ipq=0;
Ipq_xi=0;
Ipq_eta=0;
Ipq_zeta=0;

if p==q && sing
    [Ipq,Ipq_xi,Ipq_eta,Ipq_zeta] = intg_sing_SGF(k,r_cp,n1,n2,n3,3,4);
    Ipq = Ipq/(2*area); 
    % The factor of 2 here is required to get good results, but it is not
    % clear why this is correct. If removed, then the w/2 below should also be
    % changed to w, but the results are then out by 2. 
    Ipq_xi = Ipq_xi/(2*area);
    Ipq_eta = Ipq_eta/(2*area);
    Ipq_zeta = Ipq_zeta/(2*area);
    % 2A factor above required to give result in simplex coordinates - see
    % [eq 31. RWG82]. (Function intg_sing_SGF returns the LHS thereof).

else
    [w,lambda] = tri_quad(quad_pts);
    w=w/2; % Rule must be correctly normalized.
    
    %r_cp
    for nn=1:quad_pts
        r_prime = lambda(nn,1)*n1 + lambda(nn,2)*n2 + lambda(nn,3)*n3; % [eq.30,RWG82]
        R_p = norm(r_cp-r_prime); % [eq.27,RWG82]
        GF  = exp(-j*k*R_p)/R_p; % Green's function
        %w(nn)
        Ipq     = Ipq+w(nn)*GF;
        %lambda(nn,1)
        Ipq_xi  = Ipq_xi+w(nn)*lambda(nn,1)*GF;
        %lambda(nn,2)
        Ipq_eta = Ipq_eta+w(nn)*lambda(nn,2)*GF;
        %Ipq_zeta = Ipq_zeta+w(nn)*lambda(nn,3)*GF; % Code gives same
        %answers as below.
    end
    Ipq_zeta = Ipq - Ipq_xi - Ipq_eta;
end



