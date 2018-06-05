function [ifbmom] = runIFBMoMsolver(Const, zMatrices, yVectors, xVectors)
    %runIFBMoMsolver
    %   Usage:
    %       [ifbmom] = runifbMoMsolver(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       zMatrices
    %           The Z-matrices data
    %       yVectors
    %           The Yrhs-vector data
    %       xVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO)
    %
    %   Output Arguments:
    %       ifbmom
    %           Structs containing NGF-MoM solution and timing data
    %
    %   Description:
    %       Runs the IFB-MoM solution based on the Z and Y data that was read / parsed.
    %       from the FEKO *.out, *.mat, *.str and *.rhs files.
    %
    %   Note:
    %       We allow for two domains (domain A and domain B) and here we consider the MoM solution
    %       applied to each of the domains. At the moment, only a single bounce is done, which
    %       essentially corresponds to only a primary and secondary basis function that is induced
    %       on the structure. For a finite array above a ground-plane no primary MBF is generated for
    %       the ground-plane (and in turn no tertiary MBF is generated on the array).
    %
    %   References:
    %      [1] M. Carr, J. L Volakis, "Domain Decomposition by Iterative Field Bouncing"
    %      [2] Zi-Liang Liu & Chao-Fu Wang, "Efficient Iterative Method of Moments—Physical Optics
    %          Hybrid Technique for Electrically Large Objects"
    %      [3] Brand, Skrivervik, Mosig, "An Iterative Scheme Solution for the Analysis of Printed
    %          Arrays"
    %      [4] Ludick, D.J., Davidson, D.B., Jakobus, U. "Analysis of finite antenna arrays in the
    %          presence of arbitrary electromagnetic structures", ICEAA 2015
    %      [5] Ludick, D.J., Postdoc note-book 1, pp. 38 - 39. ()
    %      [6] Anastasis C. Polycarpou, "EVALUATION OF STATIONARY BLOCK ITERATIVE TECHNIQUES FOR THE
    %          SOLUTION OF FINITE ARRAYS USING THE FE-BI METHOD AND DOMAIN DECOMPOSITION"
    %      [7] Ludick, D.J, Botha, M.M., R. Maaskant, D. Davidson, "Comparison of Iterative Methods
    %          for Array Analysis,", AWPL submission 2015 (work in progress) -- Discontinued! as
    %          the DGFM-enh. Jacobi Method is not working.
    %      [8] R. Maaskant, Ludick, D.J, Botha, M.M., D. Davidson, "Analysis of Finite Antennas Arrays
    %          using the CBFM enhanced Jacobi Method,", AWPL submission 2016 (work in progress)

    error(nargchk(4,4,nargin));

    message(Const,' ');
    message(Const,'------------------------------------------------------------------------------------');
    message(Const,sprintf('Running IFB-MoM solver'));
    message(Const,sprintf('Total basis functions   : %d.',Const.numMoMbasis));
    message(Const,sprintf('Domain A basis functions: %d.',Const.numIFBbasisDomA));
    message(Const,sprintf('Domain B basis functions: %d.',Const.numIFBbasisDomB));
    message(Const,sprintf('Number of array elements: %d.',Const.numArrayElements));
    message(Const, sprintf('IFB algorithm: %d.',Const.IFBalg));
    message(Const, sprintf('IFB bounces: %d.',Const.IFB_iterations));
    message(Const, sprintf('Convergence threshold: %f percent.',Const.IFB_convergence_threshold_percentage));
    if (Const.useACA)
        message(Const, sprintf('Using ACA with tolerance %f',Const.useACAtol));
    end%if
    if ((Const.IFBalg == 12) || (Const.IFBalg == 13)|| (Const.IFBalg == 14))
        message(Const, sprintf('Number of CBFs to take into account : %d.',Const.IFB_CBFs));
        message(Const, sprintf('SVD Threshold (-1 means keep all CBFs) : %d.',Const.MBFthreshold));
        message(Const, sprintf('Use caching : %d.',Const.cache_Z0_V0));

        % Build here the file-names to reflect the particular solution configuration
        img_fcd_filename = sprintf('ifbalg%d_ifbcbfs_%d_svdthresh_%f_cacheZ0V0_%d_DGFMstart_%d_useACA_%d',...
            Const.IFBalg, Const.IFB_CBFs, Const.MBFthreshold, Const.cache_Z0_V0, Const.use_DGFM_start, Const.useACA);
    else
        img_fcd_filename = sprintf('ifbmom_ifbalg%d_useACA_%d', Const.IFBalg, Const.useACA);
    end%if

    % Set here a local debug flag on or off
    local_debug_flag = false;
    calculateRelativeResiduum = true;

    % Initialisations
    ifbmom  = [];
    ifbmom.name = 'ifbmom';
    Ntot  = Const.numMoMbasis;
    NdomA = Const.numIFBbasisDomA;
    Narr  = Ntot-NdomA;
    numSols = xVectors.numSols;             % The number of reference solutions
    ifbmom.numSols = numSols;               % Calculate a solution for each configuration
    numArrayElements = Const.numArrayElements;
    Ndgfm = Narr/numArrayElements;
    Nloc = Ndgfm;                           % Local number of basis functions for each array element
    ifbmom.Isol = complex(zeros(Ntot,1));

    % Set a flag here that is used in some cases to specify whether this is a Jacobi method 
    % (e.g. IFBalg = 12 with # CBFs = 0)
    use_Jacobi = false;

    % Set the convergence threshold here (when Const.IFB_iterations = -1)
    eps_percent = Const.IFB_convergence_threshold_percentage; 

    if (Const.IFB_iterations == -1)
        k_iter = 10000; % Set very large, essentially infinite.
    else
        k_iter = Const.IFB_iterations;
    end%if

    % Define number of CBFs to use. TO-DO: Change the name of K_back to IFB_CBFs - once M.M. Botha has seen the
    % effect of this.
    if (Const.IFB_CBFs == -1)
        K_back = k_iter;
    else
        K_back = Const.IFB_CBFs;
    end%if

    if (ifbmom.numSols > 1)
        message(Const,'runIFBMoMsolver error: Cannot run IFB-MoM solver with more than 1 sol. config.');
        error ('runIFBMoMsolver error: Cannot run IFB-MoM solver with more than 1 sol. config.');
    end

    % Start timing
    ifbmom.solTime = 0;
    tic

    % ----------------------------------------------------------------------------------------------
    if ((Const.IFBalg == 6) || (Const.IFBalg == 8))
        % For Const.IFBalg = 6 or 8 we are going to use the DGFM to calculate the starting vector
        % Setup everything so that the DGFM can be activated later when a solution is required

        % Store the current values so that we can reset them after this solver is executed
        Const.runDGFMfromIFBMoMsolver = true;
        DGFMweightVectorCalcScheme_tmp = Const.DGFMweightVectorCalcScheme;
        useDGFMmethod_tmp = Const.useDGFMmethod;
        runJACKITfromDGFM_tmp = Const.runJACKITfromDGFM;
        yVectors_tmp = yVectors;
        runNGFenDGFMsolver_tmp = Const.runNGFenDGFMsolver;

        % Change some of the configuration values (see also below for others, e.g. the weighting vector
        % scheme)
        % Only update the weighting vector scheme below (for the first iteration, i.e. the primary term,
        % we use the normal scheme=1 and therafter for each field bounce we use scheme=4)
        %Const.DGFMweightVectorCalcScheme = 4; % Use the i-DGFM with the weight-vectors the newly calculated
                                               % fields (see below where yVectors are being overridden)
        Const.useDGFMmethod = 1;               % Run on element level (not on whole sparse solution)
        Const.runJACKITfromDGFM = false;       % i-DGFM (true) or normal DGFM (false)
        Const.runNGFenDGFMsolver = false;

        % Const.useACA (not yet tested)
    end%if

    % ----------------------------------------------------------------------------------------------

    % ----------------------------------------------------------------------------------------------
    % Extract the 4 submatrices (and factorisations of the self-interaction matrices (ZdomA, ZdomB)
    % ----------------------------------------------------------------------------------------------
    % 2017-06-01: Note: We cannot apply the ACA here to calculate these submatrices. Deactivate it
    useACAtmp = Const.useACA;
    Const.useACA = 0;

    ObservRWGs = [1:NdomA];
    SourceRWGs = [1:NdomA];
    ZdomA = calcZmn(Const, zMatrices, 1, 1, ObservRWGs, SourceRWGs);

    ObservRWGs = [NdomA+1:Ntot];
    SourceRWGs = [NdomA+1:Ntot];
    ZdomB = calcZmn(Const, zMatrices, 1, 1, ObservRWGs, SourceRWGs);

    % TO-DO: Theoretically the ACA can also be used to calculate the coupling between domains A and B

    ObservRWGs = [1:NdomA];
    SourceRWGs = [NdomA+1:Ntot];
    ZdomAdomB = calcZmn(Const, zMatrices, 1, 1, ObservRWGs, SourceRWGs);
    %Zsd = zMatrices.values(1:Nngf,Nngf+1:Nmom);
    %ZdomAdomB = zMatrices.values(1:NdomA,NdomA+1:Ntot);

    ObservRWGs = [NdomA+1:Ntot];
    SourceRWGs = [1:NdomA];
    ZdomBdomA = calcZmn(Const, zMatrices, 1, 1, ObservRWGs, SourceRWGs);
    %ZdomBdomA = zMatrices.values(NdomA+1:Ntot,1:NdomA);

    Const.useACA = useACAtmp;

    % ------------------------------------
    % Calculate (ZdomA)^(-1)
    % ------------------------------------
    [LdomA,UdomA] = lu(ZdomA);

    % ------------------------------------
    % Calculate (ZdomB)^(-1)
    % ------------------------------------
    if ( (Const.IFBalg == 4) || (Const.IFBalg == 5) || (Const.IFBalg == 6) || 
         (Const.IFBalg == 7) || (Const.IFBalg == 8) || (Const.IFBalg == 9) || 
         (Const.IFBalg == 10) || (Const.IFBalg == 11) || (Const.IFBalg == 12) || ...
         (Const.IFBalg == 13) || (Const.IFBalg == 14) || (Const.IFBalg == 15))
        % Each element in domain B is considered a domain.
        % For now, we work with identical elements, so just calculate the LU factorisation
        % of the first one (also to be used for the others)
        domain_bot = Const.arrayElBasisFunctRange(1,1);
        domain_top = Const.arrayElBasisFunctRange(1,2);
        % The ACA compression of the self-term here and error in buildMoMblockACA (use therefore the full MoM submatrix)
        useACAtmp = Const.useACA;
        Const.useACA = 0;
        [ZdomB, UacadomB, VacadomB] = ...
            calcZmn(Const,zMatrices,1,1,[domain_bot:domain_top],[domain_bot:domain_top]);
        [LdomB,UdomB] = lu(ZdomB);
        % Reset again the ACA usage:
        Const.useACA = useACAtmp;
    else
        % We are considering the whole domain B as a single domain
        [LdomB,UdomB] = lu(ZdomB);
    end%if

    % ------------------------------------
    % Calculate VdomA and VdomB (initial excitation vectors)
    % ------------------------------------
    VdomA = yVectors.values(1:NdomA,1);
    VdomB = yVectors.values(NdomA+1:Ntot,1); % Note: For IFB algorithm 4 or 5, this will be overwritten

    % ----------------------------------------------------------------------------------------------
    % Start now the iterative method for calculating the total current
    % ----------------------------------------------------------------------------------------------

    % ------------------------------------
    % Iteration n_iter=0
    % ------------------------------------

    % ------------------------------------
    % Calculate JdomA_0 = (ZdomA)^(-1)*[VdomA] {The zero'th or initial iteration value - primary MBF}
    % ------------------------------------
    if (NdomA ~= 0)
        b = LdomA\VdomA;
        IdomA_0 = UdomA\b;
        ifbmom.Isol(1:NdomA,1) = IdomA_0;
    end%if

    % ------------------------------------
    % Calculate JdomB_0 = (ZdomB)^(-1)*[VdomB]
    % ------------------------------------
    if ((Const.IFBalg ~= 4) && (Const.IFBalg ~= 5) && (Const.IFBalg ~= 6) && (Const.IFBalg ~= 7) && ...
        (Const.IFBalg ~= 8) && (Const.IFBalg ~= 9) && (Const.IFBalg ~= 10) && (Const.IFBalg ~= 11)&& ...
        (Const.IFBalg ~= 12) && (Const.IFBalg ~= 13) && (Const.IFBalg ~= 14) && (Const.IFBalg ~= 15) )
        % For all algorithms other than 4/5, we can solve for whole domain B at a time
        % Alg. 4/5/6 is treated below (each array element is a domain and IdomB_0 for each is calculated
        % in a loop). The same applies to IFBalg = 5 - only then no external domain, as for IFBalg=4
        % is considdered
        b = LdomB\VdomB;
        IdomB_0 = UdomB\b;
        ifbmom.Isol(NdomA+1:Ntot,1) = IdomB_0;
    end%if

     % Stop pre-computation timing
    ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    if (Const.IFBalg == 1) % Iterative Field-bouncing algorithm as described in [1]

        % Start timing for this algorithm
        tic

        message(Const,'[runIFBMoMsolver] WARNING: IFB Algorithm 1 may lead to inaccurate results');

        % ------------------------------------
        % Iterations n_iter =1, ...
        % ------------------------------------
        % For now disable iterations (the error increases. We will need to investigate why this is the
        % case)
        for n_iter=1:1 % 1: Secondary MBFs, 2: Tertiary MBFs, etc.

            % ------------------------------------
            % Correct the minus sign for the coupling matrices (note: This was not done in [1], but
            % if the scheme is rewritten in terms of Jacobi iterations - which it actually is, then
            % the minus sign also needs to be corrected based on the particular iteration)
            fac = (-1.0d0)^n_iter;

            % ------------------------------------
            % Calculate the modified excitation vector on domain B: VdomB = - ZdomBdomA*JdomA_(n-1)
            % ------------------------------------
            VdomB = ZdomBdomA*ifbmom.Isol(1:NdomA,1);

            % ------------------------------------
            % Calculate JdomB_n = JdomB_(n-1) + (ZdomB)^(-1)*[VdomB]
            % ------------------------------------
            b = LdomB\VdomB;
            ifbmom.Isol(NdomA+1:Ntot,1) = ifbmom.Isol(NdomA+1:Ntot,1) + fac*(UdomB\b);

            % ------------------------------------
            % Calculate the modified excitation vector on domain A using the calculated current on B:
            % VdomA = - ZdomAdomB*JdomB_(n-1)
            % ------------------------------------
            VdomA = ZdomAdomB*ifbmom.Isol(NdomA+1:Ntot,1);

            % ------------------------------------
            % Calculate JdomA_n = JdomA_(n-1) + (ZdomA)^(-1)*[VdomA] {The n'th iteration value}
            % ------------------------------------
            % Back-wards substitution
            b = LdomA\VdomA;
            ifbmom.Isol(1:NdomA,1) = ifbmom.Isol(1:NdomA,1) + fac*(UdomA\b);

        end%for

         % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    elseif (Const.IFBalg == 2) % Coupling for a very specific case where we have an array (domain B)
                              % + ground plane (domain A). Only the array is excited - see ICEAA 2015
                              % article and also (postdoc small workbook 1, pp. 22-25). This then
                              % corresponds to the fact that JdomA_0 needs to be zero. The derivation
                              % is done from work by Brand, et al in [3].

        % Start timing for this algorithm
        tic

        % JdomB_0 = (ZdomB)^(-1) * VdomB        ... (1a) and
        % JdomA_0 = (ZdomA)^(-1) * VdomA = {0}  ... (1b) {Already calculated} ...

        % Make sure (1b) is satisfied (otherwise this is an error):
        if (~(all(IdomA_0(:) == 0)))
            message(Const,'[runIFBMoMsolver] IFB Algorithm cannot be used if domain A is excited');
            error ('[runIFBMoMsolver] IFB Algorithm cannot be used if domain A is excited');
        end%if

        % JdomA = [-(ZdomA)^(-1) * ZdomAdomB * JdomB_0]                                     ... (2a)
        % JdomB = [JdomB_0 + (ZdomB)^(-1) * ZdomBdomA * (ZdomA)^(-1) * ZdomAdomB * JdomB_0] ... (2b)

        % Substituting (2a) into the last term in (2b) we get:
        % JdomB = [JdomB_0 - (ZdomB)^(-1) * ZdomBdomA * JdomA]                               ... (3)

        % --------------------------------------------
        % Domain A (secondary MBF)
        % --------------------------------------------
        %IdomA_1 = inv(ZdomA)*ZdomAdomB*IdomB_0; % Rather use L,U of domain A
        b = LdomA\(ZdomAdomB*IdomB_0);
        IdomA_1 = UdomA\b;
        %ifbmom.Isol(1:NdomA,1) = -IdomA_1;  % ... (2a)

        % --------------------------------------------
        % Domain B:
        % --------------------------------------------
        b = LdomB\(ZdomBdomA*IdomA_1); % ... (primary and tertiary basis function)
        IdomB_2 = UdomB\b;
        ifbmom.Isol(NdomA+1:Ntot,1) = IdomB_0 +  IdomB_2;%  ... (3)

        % --------------------------------------------
        % Do one more correction on Domain A
        % --------------------------------------------
        b = LdomA\(ZdomAdomB*IdomB_2);
        IdomA_3 = UdomA\b;
        ifbmom.Isol(1:NdomA,1) = -IdomA_1 - IdomA_3;% ... (2a)

        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    elseif (Const.IFBalg == 3) % Essentially the same as Const.IFBalg=2, but here we support
                              % K number of iterations. If K==4, then this is exactly as
                              % Const.IFBalg=2

        % Start timing for this algorithm
        tic

        k_iter = Const.IFB_iterations; % Value > 0 (-1 not yet supported)

        Ik = complex(zeros(Ntot,k_iter));
        ifbmom.Isol(:,1) = 0.0d0; % Initialise the solution vector

        % We already have all the starting basis functions (IdomA_0 and IdomB_0 for k=1)
        Ik(1:NdomA,1) = 0.0d0;         % Domain A (primary MBF) = {0}
        Ik(NdomA+1:Ntot,1) = IdomB_0;  % Domain B (primary MBF)

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this

        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
        end%if

        for k=2:k_iter % Iterations k=2, ...
            if (mod(k,2) == 0)
                %number is even (assigned to domain A and current on domain B is zero)
                b = LdomA\(ZdomAdomB*Ik(NdomA+1:Ntot,k-1));
                Ik(1:NdomA,k) = UdomA\b;                  % Ik_domA
                Ik(NdomA+1:Ntot,k) = 0;                   % Ik_domB
                ifbmom.Isol(1:NdomA,1) = ifbmom.Isol(1:NdomA,1) - Ik(1:NdomA,k);
            else
                %number is odd (assigned to domain B and current on domain A is zero)
                Ik(1:NdomA,k) = 0;                       % Ik_domA
                b = LdomB\(ZdomBdomA*Ik(1:NdomA,k-1));
                Ik(NdomA+1:Ntot,k) = UdomB\b;            % Ik_domB
                ifbmom.Isol(NdomA+1:Ntot,1) = ifbmom.Isol(NdomA+1:Ntot,1) + Ik(NdomA+1:Ntot,k);
            end%if mod(..)
            if (Const.IFB_debug>=1)
                ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
            end%if
        end%for

        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    elseif (Const.IFBalg == 4) % See [4] and [5].

        % Start timing for this algorithm
        tic

        message(Const,'[runIFBMoMsolver] WARNING: IFB Algorithm 4 may lead to inaccurate results');

        k_iter = Const.IFB_iterations; % Value > 0 (-1 not yet supported)

        Ik = complex(zeros(Ntot,k_iter));
        ifbmom.Isol(:,1) = 0.0d0; % Initialise the solution vector

        % Domain A (primary MBF /  starting vector) .. IdomA_0 already calculated above
        Ik(1:NdomA,1) = IdomA_0;

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1)
        % for domain Bx
        % TO-DO: Danie, by improving the accuracy of the starting vector one might obtain an improved
        %        convergence rate
        for n=1:numArrayElements
            domain_bot_n = Const.arrayElBasisFunctRange(n,1);
            domain_top_n = Const.arrayElBasisFunctRange(n,2);
            VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
            b = LdomB\VdomB;
            IdomB_0 = UdomB\b;
            Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF
        end%for

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1. Note, element 1 is not part of the array, it represents the external environment,
        % e.g. a ground plane.
        if (k_iter > 1)
            for n=1:numArrayElements+1 % +1 as we might have an external domain A (e.g. ground plane)
                % Get the correct basis function range for element n:
                if (n==1)
                    domain_bot_n = 1;
                    domain_top_n = NdomA;
                    % Extract Znn^(-1) = that of external environment
                    Udom_n = UdomA;
                    Ldom_n = LdomA;
                else
                    domain_bot_n = Const.arrayElBasisFunctRange(n-1,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n-1,2);
                    % Extract Znn^(-1) = that of array element (with all
                    % elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;
                end%if
                for k=2:k_iter % Iterations k= 2, ...
                    for m=1:numArrayElements+1 % +1 as we might have an external domain A (e.g. ground plane)
                        % Add the pth coupling term for array element n using the (p-1)th coupling
                        % terms from the other array elements (discarding offcourse m=n):
                        if (m==1)
                            % The external environment
                            domain_bot_m = 1;
                            domain_top_m = NdomA;
                        else
                            % The array environment (note still access the
                            % following arrays from 1,..., numArrayElements
                            domain_bot_m = Const.arrayElBasisFunctRange(m-1,1);
                            domain_top_m = Const.arrayElBasisFunctRange(m-1,2);
                        end%if
                        if (m~=n)
                            % Extract first Znm
                            [Znm, Unm, Vnm] = ...
                                calcZmn(Const,zMatrices,1,1,[domain_bot_n:domain_top_n], ...
                                                            [domain_bot_m:domain_top_m]);
                            if (~Const.useACA)
                                b = Ldom_n\(Znm*Ik(domain_bot_m:domain_top_m,k-1));
                            else
                                % By using the ACA we can accelerate the matrix-vector product
                                % By replacing Znm ~ Unm*Vnm
                                b = Ldom_n\(Unm*Vnm*Ik(domain_bot_m:domain_top_m,k-1));
                            end%if
                            Ik(domain_bot_n:domain_top_n,k) = Udom_n\b;

                            % Add now this pth contribution to the current on element n
                            % Check again
                            ifbmom.Isol(domain_bot_n:domain_top_n,1) = ...
                                ifbmom.Isol(domain_bot_n:domain_top_n,1) + ...
                                (-1)^(k-1)*Ik(domain_bot_n:domain_top_n,k);
                        end%if
                    end%for m=1:numArrayElements

                end%for k=2:k_iter
            end%for n=1:numArrayElements
        end%if

        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    elseif ( (Const.IFBalg == 5) || (Const.IFBalg == 6)) % Similar to IFBalg = 4, but with no external domain

        % Start timing for this algorithm
        tic

        k_iter = Const.IFB_iterations; % Value > 0 (-1 not yet supported)

        Ik = complex(zeros(Ntot,k_iter));
        ifbmom.Isol(:,1) = 0.0d0; % Initialise the solution vector

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1)
        % for domain Bx
        if (Const.IFBalg == 5)
            for n=1:numArrayElements
                domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                domain_top_n = Const.arrayElBasisFunctRange(n,2);
                VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
                b = LdomB\VdomB;
                IdomB_0 = UdomB\b;
                Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF
            end%for
        else
            % Improving the accuracy of the starting vector one might obtain an improved
            % convergence rate using the DGFM
            % Use the DGFM for this step
            dgfm_sol = [];
            Const.DGFMweightVectorCalcScheme = 2; % 0 : Assume unit excitation (1+0*j)
                                                  % 1 : Standard (use VdomB)
                                                  % 2 : Use exact Xsol
            %[dgfm_sol] = runDGFMsolver(Const, zMatrices, yVectors, xVectors, [], []);
            %dgfm_sol = xVectors.values;
            %Ik(:,1) = dgfm_sol.Isol;
            Ik(:,1) = xVectors.values;

        end%if

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this
        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
        end%if

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1.
        if (k_iter > 1)
            for k=2:k_iter % Iterations k= 2, ...
                for n=1:numArrayElements
                    % Get the correct basis function range for element n:
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    % Extract Znn^(-1) = that of array element (with all elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;
                    for m=1:numArrayElements
                        if (m~=n)
                            % Add the pth coupling term for array element n using the (p-1)th coupling
                            % terms from the other array elements (discarding offcourse m=n)
                            domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                            domain_top_m = Const.arrayElBasisFunctRange(m,2);
                            % Extract first Znm
                            [Znm, Unm, Vnm] = ...
                                calcZmn(Const,zMatrices,n,m,[domain_bot_n:domain_top_n], ...
                                                            [domain_bot_m:domain_top_m]);
                            if (~Const.useACA)
                                b = Ldom_n\(Znm*Ik(domain_bot_m:domain_top_m,k-1));
                            else
                                % By using the ACA we can accelerate the matrix-vector product
                                % By replacing Znm ~ Unm*Vnm
                                b = Ldom_n\(Unm*Vnm*Ik(domain_bot_m:domain_top_m,k-1));
                            end%if
                            % I think the bug is here!
                            Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                        end%if (m~=n)
                    end%for m=1:numArrayElements

                    % Add now this pth contribution to the current on element n
                    ifbmom.Isol(domain_bot_n:domain_top_n,1) = ...
                        ifbmom.Isol(domain_bot_n:domain_top_n,1) + ...
                        (-1)^(k-1)*Ik(domain_bot_n:domain_top_n,k);
                end%for n=1:numArrayElements

                if (Const.IFB_debug>=1)
                    ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
                end%if
            end%for nk=2:k_iter

        end%if

        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    % See [3], [6] and [8] (Algs 1 and 2) ---> Rob Maaskant suggestion (not M.B. Botha -> latter move to separate block - see Const.IFBalg == 12 below.)
    elseif ((Const.IFBalg == 7)||(Const.IFBalg == 8)||(Const.IFBalg == 9)||(Const.IFBalg == 10)||(Const.IFBalg == 11))

        % Start timing for this algorithm
        tic

        Ik = complex(zeros(Ntot,k_iter));
        ifbmom.Isol(:,1) = 0.0d0; % Initialise the solution vector

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1)
        % for domain Bx
        if ((Const.IFBalg == 7)||(Const.IFBalg == 9)||(Const.IFBalg == 10)||(Const.IFBalg == 11))

            if (Const.IFBalg == 11)
                % Create first a column augmented vector of the primary
                % CBFs using the current calculated at the current iteration:
                Icbfs0 = complex(zeros(Nloc,numArrayElements)); % <-- Nloc = Ndgfm (i.e. array RWGs)
                Z0 = complex(zeros(1*numArrayElements, Ntot));  % Use only 1 CBF as testing vector
                V0 = complex(zeros(1*numArrayElements, Nloc));  % Use only 1 CBF as testing vector
            end%if

            for n=1:numArrayElements
                domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                domain_top_n = Const.arrayElBasisFunctRange(n,2);
                VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
                b = LdomB\VdomB;
                IdomB_0 = UdomB\b;
                Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF

                % For the CBFM-enh. Jacobi Method (with the precomputation step) precalculate here
                % the Z0 and V0 terms using the primary CBFs (see Alg.2 in [8])
                if (Const.IFBalg == 11)
                    Icbfs0(:,n) = IdomB_0;
                    % -- Use newly calculated CBFs for the testing vectors to precompute Z0 and V0
                    ObservRWGs = [domain_bot_n:domain_top_n];
                    SourceRWGs = [1:Ntot];
                    Z0(n,:) = ((Icbfs0(:,n)).')*calcZmn(Const, zMatrices, 1, 1, ObservRWGs, SourceRWGs);
                    V0(n,:) = ((Icbfs0(:,n)).')*yVectors.values(domain_bot_n:domain_top_n);
                end%if
            end%for

        elseif (Const.IFBalg == 8)
            % Improving the accuracy of the starting vector one might obtain an improved
            % convergence rate using the DGFM
            % Use the DGFM for this step
            dgfm_sol = [];
            Const.DGFMweightVectorCalcScheme = 1; % 0 : Assume unit excitation (1+0*j)
                                                  % 1 : Standard (use VdomB)
                                                  % 2 : Use exact Xsol
            [dgfm_sol] = runDGFMsolver(Const, zMatrices, yVectors, xVectors, [], []);

            Ik(:,1) = dgfm_sol.Isol;
        end%if

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this
        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
        end%if

        beta_k = complex(zeros(numArrayElements,1));
        beta_k(:) = 1+1i*0; % Step 1. in [7]

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1.
        if (k_iter > 1)
            for k=2:k_iter % Iterations k= 2, ...

                if (Const.IFBalg == 10)
                    % -- First store the MoM solution (will be reset below)
                    MoM_sol= xVectors.values(:,1);
                    DGFMweightVectorCalcScheme_tmp = Const.DGFMweightVectorCalcScheme;

                    % -- Overwrite now the Xsol values with that of the previous iterations
                    % -- Jacobi solution
                    xVectors.values(:,1) = Ik(:,k-1);

                    % -- Calculate now the weighting vectors, using Scheme 2 (i.e. use Xsol values)
                    Const.DGFMweightVectorCalcScheme = 2;
                    [weightVectors, tmp] = calcDGFMweightVectors(Const, zMatrices, yVectors, xVectors, []);

                    % Reset again the Xsol values (we want use them when calculating the rel. error)
                    % and also the DGFM weighting scheme
                    xVectors.values(:,1) = MoM_sol;
                    Const.DGFMweightVectorCalcScheme = DGFMweightVectorCalcScheme_tmp;
                end%if

                for n=1:numArrayElements
                    % Get the correct basis function range for element n:
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    % Extract Znn^(-1) = that of array element (with all elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;

                    for m=1:numArrayElements

                        % Add the pth coupling term for array element n using the (p-1)th coupling
                        % terms from the other array elements (discarding offcourse m=n)
                        domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                        domain_top_m = Const.arrayElBasisFunctRange(m,2);

                        % Extract first Znm
                        [Znm, Unm, Vnm] = ...
                            calcZmn(Const,zMatrices,n,m,[domain_bot_n:domain_top_n], ...
                                                        [domain_bot_m:domain_top_m]);

                        if ((Const.IFBalg == 7)||(Const.IFBalg == 8)||(Const.IFBalg == 9)||(Const.IFBalg == 11))

                            if (m~=n)
                                if (~Const.useACA)
                                    b = Ldom_n\(beta_k(m)*Znm*Ik(domain_bot_m:domain_top_m,k-1));
                                else
                                    % By using the ACA we can accelerate the matrix-vector product
                                    % By replacing Znm ~ Unm*Vnm
                                    b = Ldom_n\(beta_k(m)*Unm*Vnm*Ik(domain_bot_m:domain_top_m,k-1));
                                end%if

                                % Update surface-current at the current iteration.
                                Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                            end%if (m~=n)

                        elseif (Const.IFBalg == 10)

                            % % alpha_mn is a matrix (see FEKDDM-1.2 and also notes made on 21.06.2013)
                            alphamn = complex(zeros(Ndgfm,Ndgfm));

                            for jj = 1:Ndgfm
                                % For debugging, we can also set this equal to unity! Then we should
                                % get out exactly the Jacobi answer
                                alphamn(jj,:) = (weightVectors(:,m,1)./weightVectors(:,n,1)).';
                                %alphamn(jj,:) = 1+1i*0;
                            end%for
                            % Note: currently ACA is not supported here.

                            b = Ldom_n\((alphamn.*Znm)*Ik(domain_bot_m:domain_top_m,k-1));
                            %b = Ldom_n\(Znm*Ik(domain_bot_m:domain_top_m,k-1));
                            Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                        end%if

                    end%for m=1:numArrayElements

                    Vdom_n = yVectors.values(domain_bot_n:domain_top_n,1);
                    % See [7] (Eq. 16) : For the i-DGFM enhanced Jacobi method, we also add another
                    % term here Zpp*J_p^(k-1). Add this to Vp in
                    if (Const.IFBalg == 10)
                        Vdom_n = Vdom_n + (ZdomB*Ik(domain_bot_n:domain_top_n,k-1));
                    end%if
                    b = Ldom_n\Vdom_n;
                    Idom_n = Udom_n\b;

                    % Add now this pth contribution to the current on element n
                    ifbmom.Isol(domain_bot_n:domain_top_n,1) = Idom_n - Ik(domain_bot_n:domain_top_n,k);
                    Ik(domain_bot_n:domain_top_n,k) = ifbmom.Isol(domain_bot_n:domain_top_n,1);
                end%for n=1:numArrayElements

                % ======================================================================================================
                % Calculate the Beta_k values here. - See [8], alg's 1 and 2.
                % ======================================================================================================

                % Step 7. in [8] : Use the current iteration's values as
                % primary MBFs to calculate a new set of Beta_k values.
                if ((Const.IFBalg == 9)||(Const.IFBalg == 11))

                    % For these algorithms, we use only the latest current at the present iteration as CBF
                    numCBFsperDomain = 1;
                    cbf_start_index = k;
                    cbf_end_index = k;

                    % Create first a column augmented vector of the primary
                    % CBFs using the surface-current calculated at the current iteration (note, we only use
                    % a single CBF per domain)
                    IcbfsP = complex(zeros(Nloc,numCBFsperDomain));

                    % Calculate here space for the reduced impedance matrix
                    Zred = complex(zeros(numCBFsperDomain*numArrayElements, numCBFsperDomain*numArrayElements));
                    Vred = complex(zeros(numCBFsperDomain*numArrayElements, 1));

                    for p = 1:numArrayElements
                        % Get the correct basis function range for element p:
                        domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                        domain_top_p = Const.arrayElBasisFunctRange(p,2);
                        ind = 0;
                        for cbf_ind = cbf_start_index:cbf_end_index
                            ind = ind+1;                            
                            IcbfsP(:,ind) = Ik(domain_bot_p:domain_top_p,cbf_ind);
                        end%for

                        % Get now the CBFs for domain q:
                        IcbfsQ = complex(zeros(Nloc,numCBFsperDomain));
                        for q = 1:numArrayElements
                             % Get the correct basis function range for element q:
                            domain_bot_q = Const.arrayElBasisFunctRange(q,1);
                            domain_top_q = Const.arrayElBasisFunctRange(q,2);
                            ind = 0;
                            for cbf_ind = cbf_start_index:cbf_end_index
                                ind = ind+1;
                                IcbfsQ(:,ind) = Ik(domain_bot_q:domain_top_q,cbf_ind);
                            end%for

                            % --------------------------------------------------------------------------
                            % Calculate the correct offset for storing the entries in the reduced matrix
                            % equation (Zpq)

                            % === Range of P (rows, i.e. testing functions)
                            numMBFsP = 0;
                            for domain = 1:(p-1)
                                numMBFsP = numMBFsP + numCBFsperDomain;
                            end%for
                            PindxStart = numMBFsP+1;
                            PindxEnd   = (PindxStart - 1) + numCBFsperDomain;
                            
                            % === Range of Q (columns, i.e. basis functions)
                            numMBFsQ = 0;
                            for domain = 1:(q-1)
                                numMBFsQ = numMBFsQ + numCBFsperDomain;
                            end%for
                            QindxStart = numMBFsQ+1;
                            QindxEnd   = (QindxStart - 1) + numCBFsperDomain;

                            if (Const.IFBalg == 11)
                                % -- Use cached values where primary CBFs at k=0 was used as testing vectors
                                Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = Z0(PindxStart:PindxEnd,domain_bot_q:domain_top_q)*(IcbfsQ);
                                Vred  = V0;
                            else
                                % -- Use newly calculated CBFs for the testing vectors
                                 % Calculate the coupling matrices and Vcoupl
                                [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                                Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);
                                Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                                Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;
                            end%if

                        end%for
                        
                    end%for p = 1, numArrayElements

                    % LU-decomposition of the Z-matrix (now included here in this time)
                    % to solve now for the unknown Beta coefficients, using Zred and Vred
                    [Lred,Ured] = lu(Zred);
                    b = Lred\Vred;
                    beta_k = Ured\b;

                end%if ((Const.IFBalg == 9)||(Const.IFBalg == 11)||(Const.IFBalg == 12))

                if (Const.IFB_debug>=1)
                    ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
                end%if

                % Check whether the rel. error is <= eps_percent - if so, then exit the loop here
                if ((Const.IFB_iterations == -1) && (ifbmom.relIterError(k) <= eps_percent))
                    k_converged = k;
                    break;
                end%if

            end%for k=2:k_iter

        end%if

        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

    % =================================================================================================================
    % See [8] (Algs 3) ---> M.M. Botha suggestion where we use all the previous current values as CBFs.
    % NOTE: Here the SVD is done globally for all elements, i.e. equal number of CBFs retained for all elements.
    elseif (Const.IFBalg == 12)

        % Set here a local flag to indicate that we will be caching the values of Z0 and V0 - see [8]
        % for more information.
        cache_Z0_V0 = Const.cache_Z0_V0;

        % The above is only valid if K_back == 1 (i.e. we use only the last and latest CBF)
        if (cache_Z0_V0 && K_back ~= 1)
            message(Const, 'IFB Algorithm 12 with caching of Z0 and V0 can only be used with a single CBF');
            error('IFB Algorithm 12 with caching of Z0 and V0 can only be used with a single CBF');
        end%if

        if (calculateRelativeResiduum)
            message(Const, 'NOTE: Calculating relative residuum (execution times might be slow)');
        end%if

        % Start timing for this algorithm
        tic

        % Store all the current values at the various iterations.
        Ik = complex(zeros(Ntot,k_iter));
        ifbmom.Isol(:,1) = 0.0d0; % Initialise the solution vector

        if (cache_Z0_V0)
            % Create first a column augmented vector of the primary
            % CBFs using the current calculated at the current iteration:
            Icbfs0 = complex(zeros(Nloc,numArrayElements)); % <-- Nloc = Ndgfm (i.e. array RWGs)
            Z0 = complex(zeros(1*numArrayElements, Ntot));  % Use only 1 CBF as testing vector
            V0 = complex(zeros(1*numArrayElements, Nloc));  % Use only 1 CBF as testing vector    
        end%if

        if (Const.use_DGFM_start)
            % Use the DGFM to calculate a better guess for the initial vector (I0).
            % Improving the accuracy of the starting vector one might obtain an improved
            % convergence rate using the DGFM
            % Use the DGFM for this step
            dgfm_sol = [];
            Const.DGFMweightVectorCalcScheme = 0; % 0 : Assume unit excitation (1+0*j)
                                                  % 1 : Standard (use VdomB)
                                                  % 2 : Use exact Xsol
            [dgfm_sol] = runDGFMsolver(Const, zMatrices, yVectors, xVectors, [], []);

            Ik(:,1) = dgfm_sol.Isol;
        end%if

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1) for each
        for n=1:numArrayElements
            % Calculate the correct offsets.
            domain_bot_n = Const.arrayElBasisFunctRange(n,1);
            domain_top_n = Const.arrayElBasisFunctRange(n,2);
            if (Const.use_DGFM_start)
                IdomB_0 = Ik(domain_bot_n:domain_top_n,1);  % Domain B (array element) - primary MBF
            else
                %Standard - use the normal I0=inv(Z0)*V0 calculation
                VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
                b = LdomB\VdomB;
                IdomB_0 = UdomB\b;
                Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF
            end%if

            % For the CBFM-enh. Jacobi Method (with the precomputation step) precalculate here
            % the Z0 and V0 terms using the primary CBFs (see Alg.2 in [8])
            if (cache_Z0_V0)
                Icbfs0(:,n) = IdomB_0;
                % -- Use newly calculated CBFs for the testing vectors to precompute Z0 and V0
                ObservRWGs = [domain_bot_n:domain_top_n];
                SourceRWGs = [1:Ntot];
                Z0(n,:) = ((Icbfs0(:,n)).')*calcZmn(Const, zMatrices, 1, 1, ObservRWGs, SourceRWGs);
                V0(n,:) = ((Icbfs0(:,n)).')*yVectors.values(domain_bot_n:domain_top_n);
            end%if
        end%for

        % Store the 0th iteration.
        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this

        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

        % Note: Zero the timing here, as all the methods start off with the same solution (equal to the
        % isolated case).
        ifbmom.solTime = 0; % Zero the value here.

        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.iterTiming = zeros(1,k_iter);
                        % Added now the relative residuum norm here also (only if enabled, as this is costly)
            if (calculateRelativeResiduum)
                ifbmom.relResError = zeros(1,k_iter);
            end
            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
            if (calculateRelativeResiduum)
                ifbmom.relResError(1) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
            end%if
            ifbmom.iterTiming(1) = ifbmom.solTime;
        end%if

        % Use the Jacobi method when K_back = 0
        use_Jacobi = false;
        if (K_back == 0)
            use_Jacobi = true;
        end%if

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1.
        if (k_iter > 1)
            for k=2:k_iter % Iterations k= 2, ...

                % Time now each iteration
                tic

                message(Const,sprintf('\nIteration k = %d',k));

                % =====================================================================================================
                % Calculate the Beta_k values here. - See [8], Alg. 3 - other than Alg. 1 and 2, the Beta coefficients
                % are calculated beforehand here.
                % =====================================================================================================
                if (~use_Jacobi)
                    if (false)
                        % -- Version 1 : Use 1 2 3 .... N
                        %                        {------} these CBFs
                        %                         K_back

                        % For this algorithm, we use only the latest current at the present iteration as CBF
                        if (k <= K_back)
                            cbf_start_index = 1;
                        else
                            cbf_start_index = k - K_back; % + 1; --> not correct to (+1).
                        end
                        cbf_end_index = k-1;

                    else
                        if (K_back == 1)
                            cbf_start_index = k-1;
                            cbf_end_index = k-1;
                        else
                            % -- Version 2 : Use 1 2 3 .... N
                            %                   {------} these CBFs ... Better, as the initial ones are more orthogonal
                            %                     K_back

                            % For this algorithm, we always start from CBF index 1
                            cbf_start_index = 1;
                            if (k <= K_back)
                                cbf_end_index = k-1;
                            else
                                cbf_end_index = K_back; % + 1; --> not correct to (+1).
                            end
                        end % if (K_back == 1)
                    end%if

                    % Update now the number of CBFs that we will use.
                    numCBFsperDomain = cbf_end_index - cbf_start_index + 1;

                    % Check now whether we apply an orthonormalisation scheme to improve the quality of the CBFs, and
                    % retain only those above a certain threshold (when applying the SVD).
                    % 2016-06-14 (originally started the 2016-01-14): 
                    %     Reduce and orthonormalize the CBFs here (see also runMBFsolver for 
                    %             similar approach).
                    if (Const.useMBFreduction && ~(K_back==1))
                        %tic --> Need to time this ... NB!!!
                        message(Const,sprintf('  Reduce and orthonormalise CBFs'));
                        
                        % Put all the (K_back #) CBFs in a column augmented matrix for SVD reduction.
                        % Note: we reduce the "global" CBfs - i.e. not per array element.
                        origCBFs = complex(zeros(Ntot,cbf_end_index));
                        message(Const,sprintf('    Number of initially generated CBFs  = %d',numCBFsperDomain));
                        ind = 0;
                        for cbf_ind = cbf_start_index:cbf_end_index
                            ind = ind+1;                            
                            origCBFs(:,ind) = Ik(:,cbf_ind);
                        end%for

                        % Plot here the SVD spectrum (only for the last, and specified iteration)
                        % and only if in local debug mode.
                        k_plot_SVD = k_iter+1;
                        plot_SVD = true;
                        if ( plot_SVD && ((k == k_iter) || (k == k_plot_SVD)) )
                            % Plot the data.
                            MBFplotSVspectrum_tmp = Const.MBFplotSVspectrum;
                            Const.MBFplotSVspectrum = true;

                            % Make sure the data is store to an *.FCD file for further processing in POSTFEKO.
                            store_to_fcd_file_tmp = Const.store_to_fcd_file;
                            Const.store_to_fcd_file = true;
                        end%if

                        % The singular values are stored to an FCD string - 
                        fcdString = sprintf('ifb_svd_alg%d',Const.IFBalg);
                        redCBFs = reduceMBFset(Const, origCBFs, fcdString);

                        % Reset the SVD spectrum plotting and FCD file storage if enabled above.
                        if (plot_SVD && ((k == k_iter) || (k == k_plot_SVD)) )
                            Const.MBFplotSVspectrum = MBFplotSVspectrum_tmp;
                            Const.store_to_fcd_file = store_to_fcd_file_tmp;
                        end%if

                        % Retain now the orthonormalised CBFs
                        Ik_ortho = redCBFs;
                        
                        % Update the number of CBFs that will be taken into account - note the _ortho arrays and vectors
                        % always run from cbf_start_index = 1!
                        numCBFsperDomain = size(redCBFs,2);
                        cbf_end_index = numCBFsperDomain;
                        cbf_start_index = 1;
                        
                        message(Const,sprintf('    Number of retained orthonormal CBFs  = %d',size(redCBFs,2)));
                                        
                        %ifbmom.svdTime(solNum) = toc; % End timing
                    else
                        %ifbmom.svdTime(solNum) = 0.0;
                    end %if (Const.useMBFreduction)

                    % Calculate here space for the reduced impedance matrix
                    Zred = complex(zeros(numCBFsperDomain*numArrayElements, numCBFsperDomain*numArrayElements));
                    Vred = complex(zeros(numCBFsperDomain*numArrayElements, 1));

                    % Create first a column augmented vector of the primary
                    % CBFs using the surface-current calculated at the current iteration(s).
                    IcbfsP = complex(zeros(Nloc,numCBFsperDomain));
                    for p = 1:numArrayElements
                        % Get the correct basis function range for element p:
                        domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                        domain_top_p = Const.arrayElBasisFunctRange(p,2);
                        ind = 0;
                        for cbf_ind = cbf_start_index:cbf_end_index
                            ind = ind+1;                            
                            if (Const.useMBFreduction && ~(K_back==1))
                                % -- Use post-SVD reduced CBFs
                                IcbfsP(:,ind) = Ik_ortho(domain_bot_p:domain_top_p,ind);
                            else
                                % -- Use standard CBFs (not orthonormalised)
                                IcbfsP(:,ind) = Ik(domain_bot_p:domain_top_p,cbf_ind);
                            end%if
                        end%for

                        % Get now the CBFs for domain q:
                        IcbfsQ = complex(zeros(Nloc,numCBFsperDomain));
                        for q = 1:numArrayElements
                             % Get the correct basis function range for element q:
                            domain_bot_q = Const.arrayElBasisFunctRange(q,1);
                            domain_top_q = Const.arrayElBasisFunctRange(q,2);
                            ind = 0;
                            for cbf_ind = cbf_start_index:cbf_end_index
                                ind = ind+1;
                                if (Const.useMBFreduction&& ~(K_back==1))
                                    % -- Use post-SVD reduced CBFs
                                    IcbfsQ(:,ind) = Ik_ortho(domain_bot_q:domain_top_q,ind);
                                else
                                    % -- Use standard CBFs (not orthonormalised)
                                    IcbfsQ(:,ind) = Ik(domain_bot_q:domain_top_q,cbf_ind);
                                end%if
                            end%for

                            % --------------------------------------------------------------------------
                            % Calculate the correct offset for storing the entries in the reduced matrix
                            % equation (Zpq)

                            % === Range of P (rows, i.e. testing functions)
                            numMBFsP = 0;
                            for domain = 1:(p-1)
                                numMBFsP = numMBFsP + numCBFsperDomain;
                            end%for
                            PindxStart = numMBFsP+1;
                            PindxEnd   = (PindxStart - 1) + numCBFsperDomain;
                            
                            % === Range of Q (columns, i.e. basis functions)
                            numMBFsQ = 0;
                            for domain = 1:(q-1)
                                numMBFsQ = numMBFsQ + numCBFsperDomain;
                            end%for
                            QindxStart = numMBFsQ+1;
                            QindxEnd   = (QindxStart - 1) + numCBFsperDomain;

                            if (cache_Z0_V0)
                                % -- Use cached values where primary CBFs at k=0 was used as testing vectors
                                Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = Z0(PindxStart:PindxEnd,domain_bot_q:domain_top_q)*(IcbfsQ);
                                Vred  = V0;
                            else
                                % -- Use newly calculated CBFs for the testing vectors
                                 % Calculate the coupling matrices and Vcoupl
                                [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                                Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);
                                Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                                Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;
                            end%if

                            % % -- Use newly calculated CBFs for the testing vectors
                            % %    Calculate the coupling matrices, Zcoupl (or Ucoupl and Vcoupl if the 
                            % %    ACA is used)
                            % [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                            % Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);
                            
                            % Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                            % Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;

                        end%for q = 1:numArrayElements
                        
                    end%for p = 1, numArrayElements

                    % LU-decomposition of the Z-matrix (now included here in this time)
                    % to solve now for the unknown Beta coefficients, using Zred and Vred
                    [Lred,Ured] = lu(Zred);
                    b = Lred\Vred;
                    beta_k = Ured\b;  % A one-dimensional array of length : [numCBFsperDomain * numArrayElements x 1]

                    if (local_debug_flag)
                        % Some debug output here.
                        fprintf('*** debug : beta_k = %f + 1i*%f\n',real(beta_k), imag(beta_k));
                    end%if

                end% if (~use_Jacobi)

                % Now that we have the beta_k coefficients, we can apply them to the solution to calculate the surface
                % current at the present iteration.

                % Switch off ACA calculation here for the mat-vec below
                useACAtmp = Const.useACA;
                Const.useACA = 0;

                for n=1:numArrayElements % To calculate the current on domain p - as noted in [8] (used "n" here just 
                                         % for consistency witht the above code)
                    % Get the correct basis function range for element n:
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    % Extract Znn^(-1) = that of array element (with all elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;

                    for m=1:numArrayElements

                        % Add the pth coupling term for array element n using the (p-1)th coupling
                        % terms from the other array elements (discarding offcourse m=n)
                        domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                        domain_top_m = Const.arrayElBasisFunctRange(m,2);

                        % Extract first Znm
                        [Znm, Unm, Vnm] = ...
                            calcZmn(Const,zMatrices,n,m,[domain_bot_n:domain_top_n], ...
                                                        [domain_bot_m:domain_top_m]);

                        if (m~=n)
                            if (~use_Jacobi)
                                ind = 0;
                                for cbf_ind = cbf_start_index:cbf_end_index
                                    ind = ind+1;
                                    % Distinguish whether we are using the orthonormalised CBFs or not
                                    if (Const.useMBFreduction && ~(K_back==1))
                                        Ik_tmp = Ik_ortho(domain_bot_m:domain_top_m,ind);
                                    else
                                        Ik_tmp = Ik(domain_bot_m:domain_top_m,cbf_ind);
                                    end%if

                                    % beta_k is stored in one-dimensional indexing. Calculate here this index, based on 
                                    % the following general formula for a 2-D array stored in 1-D formatting:
                                    %       A(I,J) = I + (J-1)*LDAMN - where LDAMN is the leading dimension. 
                                    % In our case, the followng values hold:
                                    %       I : ind (i.e. the basis function index)
                                    %       J : m (i.e. the element number)
                                    %       LDAMB : numCBFsperDomain (number of CBFs / domain - here equal for all 
                                    %               array elements. See Alg. 13 where this might not be the case.)
                                    beta_k_1d_index = ind+(m-1)*numCBFsperDomain;
                                    if (~Const.useACA)
                                        if (false)
                                            % Some debug output here.
                                            fprintf('*** debug : size(Znm) = [%d,%d]\n',size(Znm,1),size(Znm,2));
                                            fprintf('*** debug : size(Ik_tmp) = [%d,%d]\n',size(Ik_tmp,1),size(Ik_tmp,2));
                                            fprintf('*** debug : size(beta_k) = [%d,%d]\n',size(beta_k,1),size(beta_k,2));
                                        end%if
                                        b = Ldom_n\(beta_k(beta_k_1d_index,1)*Znm*Ik_tmp);
                                    else
                                        % By using the ACA we can accelerate the matrix-vector product
                                        % By replacing Znm ~ Unm*Vnm
                                        b = Ldom_n\(beta_k(beta_k_1d_index,1)*Unm*Vnm*Ik_tmp);
                                    end%if

                                    % Update surface-current at the current iteration.
                                    Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                                end%for cbf_ind = cbf_start_index:cbf_end_index
                            else
                                % We do are using the Jacobi method - i.e. no Beta coefficients used.
                                if (~Const.useACA)
                                    b = Ldom_n\(Znm*Ik(domain_bot_m:domain_top_m,k-1));
                                else
                                    % By using the ACA we can accelerate the matrix-vector product
                                    % By replacing Znm ~ Unm*Vnm
                                    b = Ldom_n\(Unm*Vnm*Ik(domain_bot_m:domain_top_m,k-1));
                                end%if

                                % Update surface-current at the current iteration.
                                Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                            end% if ~(use_Jacobi)

                        end%if (m~=n)
                    end%for m=1:numArrayElements

                    Vdom_n = yVectors.values(domain_bot_n:domain_top_n,1);
                    b = Ldom_n\Vdom_n;
                    Idom_n = Udom_n\b;

                    % Add now this pth contribution to the current on element n
                    ifbmom.Isol(domain_bot_n:domain_top_n,1) = Idom_n - Ik(domain_bot_n:domain_top_n,k);
                    Ik(domain_bot_n:domain_top_n,k) = ifbmom.Isol(domain_bot_n:domain_top_n,1);
                end%for n=1:numArrayElements

                % See above comment - switch the ACA back to its original status
                Const.useACA = useACAtmp;  

                % Stop pre-computation timing
                ifbmom.solTime = ifbmom.solTime + toc;

                if (Const.IFB_debug>=1)
                    ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
                    % Calculate the relative residual error here (only if enabled, as this is expensive)
                    if (calculateRelativeResiduum)
                        ifbmom.relResError(k) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
                    end%if
                    ifbmom.iterTiming(k) = ifbmom.solTime;
                    if (local_debug_flag)
                        % Some debug output here.
                        fprintf('*** debug : ifbmom.relIterError(k) = %f %%\n',ifbmom.relIterError(k));
                    end%if
                end%if

                % Check whether the rel. error is <= eps_percent - if so, then exit the loop here
                if ((Const.IFB_iterations == -1) && (ifbmom.relIterError(k) <= eps_percent))
                    k_converged = k;
                    break;
                end%if

            end%for k=2:k_iter

        end%if

    % =================================================================================================================
    % See [8] (Algs 3) ---> M.M. Botha suggestion where we use all the previous current values as CBFs.
    % NOTE: Here the SVD is done per element, i.e. the number of CBFs retained for each element might be different.
    elseif (Const.IFBalg == 13)

        if (~Const.useMBFreduction)
            message(Const, 'IFB Algorithm 13 can only be used with SVD reduction');
            error('IFB Algorithm 13 can only be used with SVD reduction');
        end%if

        if (Const.IFB_CBFs==1)
            message(Const, 'IFB Algorithm 13 can only be used when all previous iterations are used as CBFs');
            error('IFB Algorithm 13 can only be used when all previous iterations are used as CBFs');
        end%if

        if (calculateRelativeResiduum)
            message(Const, 'NOTE: Calculating relative residuum (execution times might be slow)');
        end%if    

        % Start timing for this algorithm
        tic

        % Store all the current values at the various iterations.
        Ik = complex(zeros(Ntot,k_iter));

        % Make sure we have enough space to store the orthonormalised eigencurrents for each array element.
        % Use a 3D matrix for this (Nrwg, # CBFs, # array el. id); Use also a vector to indicate how many orthogonal
        % eigencurrents there are per array element.
        Ik_ortho = complex(zeros(Ntot,K_back,numArrayElements));
        num_orth_CBFs_per_element = zeros(numArrayElements,1);
        
        % Initialise the solution vector
        ifbmom.Isol(:,1) = 0.0d0;

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1) for each
        for n=1:numArrayElements
            domain_bot_n = Const.arrayElBasisFunctRange(n,1);
            domain_top_n = Const.arrayElBasisFunctRange(n,2);
            VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
            b = LdomB\VdomB;
            IdomB_0 = UdomB\b;
            Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF
        end%for

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this
        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

        % Note: Zero the timing here, as all the methods start off with the same solution (equal to the
        % isolated case).
        ifbmom.solTime = 0; % Zero the value here.
                                       
        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.iterTiming = zeros(1,k_iter);

            % Also profile now the reduced matrix setup and the orthonormalisation step
            ifbmom.orthonormTiming = zeros(1,k_iter);
            ifbmom.cbfmTiming = zeros(1,k_iter);
            
            % Added now the relative residuum norm here also (only if enabled, as this is costly)
            if (calculateRelativeResiduum)
                ifbmom.relResError = zeros(1,k_iter);
            end
            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
            if (calculateRelativeResiduum)
                ifbmom.relResError(1) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
            end%if
            ifbmom.iterTiming(1) = ifbmom.solTime;
            ifbmom.orthonormTiming(1) = 0;
            ifbmom.cbfmTiming(1) = 0;
        end%if

        %beta_k = complex(zeros(numArrayElements,K_back)); % We actually store only K (Zeta) number of Beta's.

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1.
        if (k_iter > 1)
            for k=2:k_iter % Iterations k= 2, ...

                % Time now each iteration
                %tic

                message(Const,sprintf('\nIteration k = %d',k));

                % =====================================================================================================
                % Calculate the Beta_k values here. - See [8], Alg. 3 - other than Alg. 1 and 2, the Beta coefficients
                % are calculated beforehand here.
                % =====================================================================================================

                if (false)
                    % -- Version 1 : Use 1 2 3 .... N
                    %                        {------} these CBFs
                    %                         K_back

                    % For this algorithm, we use only the latest current at the present iteration as CBF
                    if (k <= K_back)
                        cbf_start_index = 1;
                    else
                        cbf_start_index = k - K_back; % + 1; --> not correct to (+1).
                    end
                    cbf_end_index = k-1;

                else
                    % -- Version 1 : Use 1 2 3 .... N
                    %                   {------} these CBFs ... Better, as the initial ones are more orthogonal
                    %                     K_back

                    % For this algorithm, we always start from CBF index 1
                    cbf_start_index = 1;
                    if (k <= K_back)
                        cbf_end_index = k-1;
                    else
                        cbf_end_index = K_back; % + 1; --> not correct to (+1).
                    end
                end%if

                % Update now the number of CBFs that we will use.
                numCBFsperDomain = cbf_end_index - cbf_start_index + 1;

                % Start timing the orthonormalisation
                tic

                % Apply now an orthonormalisation scheme to improve the quality of the CBFs, and retain only those 
                % above a certain threshold (when applying the SVD). (See also runMBFsolver for similar approach).
                
                % Loop over all the elements and reduce the CBF set of each individually (other than Alg. 12 where
                % this is done globally).
                totRedCBFs = 0;
                % Put all the (K_back #) CBFs in a column augmented matrix for SVD reduction.
                % Note: we reduce the "global" CBfs - i.e. not per array element.
                origCBFs = complex(zeros(Nloc,numCBFsperDomain));
                for p = 1:numArrayElements
                    % Get the correct basis function range for element p:
                    domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                    domain_top_p = Const.arrayElBasisFunctRange(p,2);

                    message(Const,sprintf('    Reduce and orthonormalise CBFs for array element %d',p));
                    
                    % TO-DO: 2016-11-13: Danie, rechecking this part now again to ensure we orthonormalise the 
                    % CBFs correctly

                    origCBFs(:,:) = 0+1i*0;
                    if (true)
                        fprintf('*** debug : Number of initially generated CBFs = %d\n',numCBFsperDomain);
                    end%if
                    ind = 0;
                    for cbf_ind = cbf_start_index:cbf_end_index
                        ind = ind+1;                            
                        origCBFs(:,ind) = Ik(domain_bot_p:domain_top_p,cbf_ind);
                    end%for

                    % Plot here the SVD spectrum (only for the last, and specified iteration)
                    % and only if in local debug mode.
                    k_plot_SVD = k_iter+1;
                    plot_SVD = false;
                    if ( plot_SVD && ((k == k_iter) || (k == k_plot_SVD)) )
                        % Plot the data.
                        MBFplotSVspectrum_tmp = Const.MBFplotSVspectrum;
                        Const.MBFplotSVspectrum = true;

                        % Make sure the data is store to an *.FCD file for further processing in POSTFEKO.
                        store_to_fcd_file_tmp = Const.store_to_fcd_file;
                        Const.store_to_fcd_file = true;
                    end%if

                    % The singular values are stored to an FCD string - 
                    fcdString = sprintf('ifb_svd_alg%d_arryEl%d',Const.IFBalg,p);
                    redCBFs = reduceMBFset(Const,origCBFs,fcdString);

                    % Reset the SVD spectrum plotting and FCD file storage if enabled above.
                    if (plot_SVD && ((k == k_iter) || (k == k_plot_SVD)) )
                        Const.MBFplotSVspectrum = MBFplotSVspectrum_tmp;
                        Const.store_to_fcd_file = store_to_fcd_file_tmp;
                    end%if

                    % Increment the total number of CBFs (all domains), so that we can allocate space for Zred and
                    % Vred correctly below.
                    totRedCBFs = totRedCBFs + size(redCBFs,2);
                    
                    % Retain now the orthonormalised CBFs (also update the number of orth. CBFs / element).
                    Ik_ortho(domain_bot_p:domain_top_p,1:size(redCBFs,2),p) = redCBFs;
                    num_orth_CBFs_per_element(p) = size(redCBFs,2);
                    
                    if (true)
                        fprintf('*** debug : Number of retained orthonormal CBFs = %d\n\n',size(redCBFs,2));
                    end%if                
                    %ifbmom.svdTime(solNum) = toc; % End timing
                end% for p = 1:numArrayElements

                % Stop the orthonormalisation timing
                ifbmom.orthonormTiming(k) = ifbmom.orthonormTiming(k-1) + toc;

                % Start timing the CBFM solution
                tic

                % Calculate here space for the reduced impedance matrix
                Zred = complex(zeros(totRedCBFs, totRedCBFs));
                Vred = complex(zeros(totRedCBFs, 1));

                % Fill now Zred and Vred to calculate the Beta_k's, using the orthonormalised CBFs
                for p = 1:numArrayElements

                    % Copy only the reduced CBFs
                    IcbfsP = complex(zeros(Nloc,num_orth_CBFs_per_element(p)));

                    % Get the correct basis function range for element p:
                    domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                    domain_top_p = Const.arrayElBasisFunctRange(p,2);
                    
                    for cbf_ind = 1:num_orth_CBFs_per_element(p)
                        % -- Use post-SVD reduced CBFs
                        IcbfsP(:,cbf_ind) = Ik_ortho(domain_bot_p:domain_top_p,cbf_ind,p);
                    end%for

                    % Get now the CBFs for domain q:
                    for q = 1:numArrayElements
                        % Copy only the reduced CBFs
                        IcbfsQ = complex(zeros(Nloc,num_orth_CBFs_per_element(q)));

                         % Get the correct basis function range for element q:
                        domain_bot_q = Const.arrayElBasisFunctRange(q,1);
                        domain_top_q = Const.arrayElBasisFunctRange(q,2);
                        for cbf_ind = 1:num_orth_CBFs_per_element(q)
                             % -- Use post-SVD reduced CBFs
                            IcbfsQ(:,cbf_ind) = Ik_ortho(domain_bot_q:domain_top_q,cbf_ind,q);
                        end%for

                        % --------------------------------------------------------------------------
                        % Calculate the correct offset for storing the entries in the reduced matrix
                        % equation (Zpq)

                        % === Range of P (rows, i.e. testing functions)
                        numMBFsP = 0;
                        for domain = 1:(p-1)
                            numMBFsP = numMBFsP + num_orth_CBFs_per_element(domain);
                        end%for
                        PindxStart = numMBFsP+1;
                        PindxEnd   = (PindxStart - 1) + num_orth_CBFs_per_element(p);
                        
                        % === Range of Q (columns, i.e. basis functions)
                        numMBFsQ = 0;
                        for domain = 1:(q-1)
                            numMBFsQ = numMBFsQ + num_orth_CBFs_per_element(domain);
                        end%for
                        QindxStart = numMBFsQ+1;
                        QindxEnd   = (QindxStart - 1) + num_orth_CBFs_per_element(q);

                        % -- Use newly calculated CBFs for the testing vectors
                        %    Calculate the coupling matrices, Zcoupl (or Ucoupl and Vcoupl if the 
                        %    ACA is used) and Vcoupl
                        [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                        Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);
                        
                        Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                        Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;

                    end%for q = 1:numArrayElements
                end%for p = 1, numArrayElements

                % calculate the memory usage here of Zred:
                ifbmom.memUsageZred = byteSize(Zred);

                % LU-decomposition of the Z-matrix (now included here in this time)
                % to solve now for the unknown Beta coefficients, using Zred and Vred
                [Lred,Ured] = lu(Zred);
                b = Lred\Vred;
                beta_k = Ured\b;  % A one-dimensional array of length : [numCBFsperDomain * numArrayElements x 1]

                if (false)
                    % Some debug output here.
                    beta_k
                    %fprintf('*** debug : beta_k = %f + 1i*%f\n',real(beta_k), imag(beta_k));
                end%if

                % Now that we have the beta_k coefficients, we can apply them to the solution to calculate the surface
                % current at the present iteration.

                for n=1:numArrayElements % To calculate the current on domain p - as noted in [8] (used "n" here just 
                                         % for consistency witht the above code)
                    % Get the correct basis function range for element n:
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    % Extract Znn^(-1) = that of array element (with all elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;

                    for m=1:numArrayElements

                        % Add the pth coupling term for array element n using the (p-1)th coupling
                        % terms from the other array elements (discarding offcourse m=n)
                        domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                        domain_top_m = Const.arrayElBasisFunctRange(m,2);

                        % Extract first Znm
                        [Znm, Unm, Vnm] = ...
                            calcZmn(Const,zMatrices,n,m,[domain_bot_n:domain_top_n], ...
                                                        [domain_bot_m:domain_top_m]);

                        if (m~=n)

                            % Calculate the correct beta_k_offset for the array element, i.e. where its beta coefficients
                            % are stored in the beta_k vector (remember, some elements might have more than others).
                            % This value is calculated by considering the number of beta coefficients of the previous
                            % array elements.
                            beta_k_offset = 0;
                            for array_el = 1:(m-1)
                                beta_k_offset = beta_k_offset + num_orth_CBFs_per_element(array_el);
                            end%for

                            for cbf_ind = 1:num_orth_CBFs_per_element(m)
                                Ik_tmp = Ik_ortho(domain_bot_m:domain_top_m,cbf_ind,m);

                                % Calculate now the correct beta index (see also comment above)
                                beta_k_index = beta_k_offset + cbf_ind;

                                if (~Const.useACA)
                                    if (false)
                                        % Some debug output here.
                                        fprintf('\n*** debug : size(Z%d_%d) = [%d,%d]\n',n,m,size(Znm,1),size(Znm,2));
                                        fprintf('*** debug : size(Ik_tmp)  = [%d,%d]\n',size(Ik_tmp,1),size(Ik_tmp,2));
                                        fprintf('*** debug : size(beta_k)  = [%d,%d]\n',size(beta_k,1),size(beta_k,2));
                                        fprintf('*** debug : beta_k_offset = %d\n',beta_k_offset);
                                        fprintf('*** debug : beta_k_index  = %d\n',beta_k_index);
                                    end%if
                                    b = Ldom_n\(beta_k(beta_k_index,1)*Znm*Ik_tmp);
                                else
                                    % By using the ACA we can accelerate the matrix-vector product
                                    % By replacing Znm ~ Unm*Vnm
                                    b = Ldom_n\(beta_k(beta_k_index,1)*Unm*Vnm*Ik_tmp);
                                end%if

                                % Update surface-current at the current iteration.
                                Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                            end%for cbf_ind = cbf_start_index:cbf_end_index
                        end%if (m~=n)
                    end%for m=1:numArrayElements

                    Vdom_n = yVectors.values(domain_bot_n:domain_top_n,1);
                    b = Ldom_n\Vdom_n;
                    Idom_n = Udom_n\b;

                    % Add now this pth contribution to the current on element n
                    ifbmom.Isol(domain_bot_n:domain_top_n,1) = Idom_n - Ik(domain_bot_n:domain_top_n,k);
                    Ik(domain_bot_n:domain_top_n,k) = ifbmom.Isol(domain_bot_n:domain_top_n,1);
                end%for n=1:numArrayElements

                % Stop the CBFM timing
                ifbmom.cbfmTiming(k) = ifbmom.cbfmTiming(k-1) + toc;


                if (Const.IFB_debug>=1)
                    ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
                    % Calculate the relative residual error here (only if enabled, as this is expensive)
                    if (calculateRelativeResiduum)
                        ifbmom.relResError(k) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
                    end%if
                    % Calculate the total time - sum of the Orthonormalisation and CBFM timing
                    ifbmom.iterTiming(k) = ifbmom.orthonormTiming(k) + ifbmom.cbfmTiming(k);
                    if (local_debug_flag)
                        % Some debug output here.
                        fprintf('*** debug : ifbmom.relIterError(k) = %f %%\n',ifbmom.relIterError(k));
                    end%if
                end%if

                % Check whether the rel. error is <= eps_percent - if so, then exit the loop here
                if ((Const.IFB_iterations == -1) && (ifbmom.relIterError(k) <= eps_percent))
                    k_converged = k;
                    break;
                end%if

            end%for k=2:k_iter

        end%if

        % Make sure we record the total time again:
        ifbmom.solTime = ifbmom.iterTiming(k_iter);

    % =================================================================================================================
    % 2016-11-13: This is the same as Alg. 12, but with the Orthonormalisation done per element (and not globally)
    elseif (Const.IFBalg == 14)
    
        if (Const.useMBFreduction)
            message(Const, 'IFB Algorithm 14 cannot be used with SVD reduction');
            error('IFB Algorithm 14 cannot be used with SVD reduction');
        end%if

        if (Const.IFB_CBFs==1)
            message(Const, 'IFB Algorithm 14 can only be used when all previous iterations are used as CBFs');
            error('IFB Algorithm 14 can only be used when all previous iterations are used as CBFs');
        end%if

        if (Const.IFB_CBFs==1)
            message(Const, 'IFB Algorithm 14 cannot be used for a Jacobi only simulation');
            error('IFB Algorithm 14 cannot be used for a Jacobi only simulation');
        end%if

        % Start timing for this algorithm
        tic

        % Store all the current values at the various iterations.
        Ik = complex(zeros(Ntot,k_iter));

        % Make sure we have enough space to store the orthonormalised eigencurrents for each array element.
        %Ik_ortho = complex(zeros(Ntot,K_back));
        Ik_ortho = complex(zeros(Ntot,k_iter));
        
        % Initialise the solution vector
        ifbmom.Isol(:,1) = 0.0d0;

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1) for each
        for n=1:numArrayElements
            domain_bot_n = Const.arrayElBasisFunctRange(n,1);
            domain_top_n = Const.arrayElBasisFunctRange(n,2);
            VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
            b = LdomB\VdomB;
            IdomB_0 = UdomB\b;
            Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF
        end%for

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this
        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

        % Note: Zero the timing here, as all the methods start off with the same solution (equal to the
        % isolated case).
        ifbmom.solTime = 0; % Zero the value here.
                                       
        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.iterTiming = zeros(1,k_iter);

            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
            ifbmom.iterTiming(1) = ifbmom.solTime;
        end%if

        %beta_k = complex(zeros(numArrayElements,K_back)); % We actually store only K (Zeta) number of Beta's.

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1.
        if (k_iter > 1)
            for k=2:k_iter % Iterations k= 2, ...

                % Time now each iteration
                tic

                message(Const,sprintf('\nIteration k = %d',k));

                % =====================================================================================================
                % Calculate the Beta_k values here. - See [8], Alg. 3 - other than Alg. 1 and 2, the Beta coefficients
                % are calculated beforehand here.
                % =====================================================================================================

                if (false)
                    % -- Version 1 : Use 1 2 3 .... N
                    %                        {------} these CBFs
                    %                         K_back

                    % For this algorithm, we use only the latest current at the present iteration as CBF
                    if (k <= K_back)
                        cbf_start_index = 1;
                    else
                        cbf_start_index = k - K_back; % + 1; --> not correct to (+1).
                    end
                    cbf_end_index = k-1;

                else
                    % -- Version 1 : Use 1 2 3 .... N
                    %                   {------} these CBFs ... Better, as the initial ones are more orthogonal
                    %                     K_back

                    % For this algorithm, we always start from CBF index 1
                    cbf_start_index = 1;
                    if (k <= K_back)
                        cbf_end_index = k-1;
                    else
                        cbf_end_index = K_back; % + 1; --> not correct to (+1).
                    end
                end%if

                % Update now the number of CBFs that we will use.
                numCBFsperDomain = cbf_end_index - cbf_start_index + 1;

                % Apply now an orthonormalisation scheme to improve the quality of the CBFs, and retain only those 
                % above a certain threshold (when applying the SVD). (See also runMBFsolver for similar approach).
                
                % Loop over all the elements and reduce the CBF set of each individually (other than Alg. 12 where
                % this is done globally).
                totRedCBFs = 0;
                % Put all the (K_back #) CBFs in a column augmented matrix for SVD reduction.
                % Note: we reduce the "global" CBfs - i.e. not per array element.
                origCBFs = complex(zeros(Nloc,numCBFsperDomain));
                for p = 1:numArrayElements
                    % Get the correct basis function range for element p:
                    domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                    domain_top_p = Const.arrayElBasisFunctRange(p,2);

                    message(Const,sprintf('    Reduce and orthonormalise CBFs for array element %d',p));
                    
                    % TO-DO: 2016-11-13: Danie, rechecking this part now again to ensure we orthonormalise the 
                    % CBFs correctly - i.e. we only orthonormalise the CBFs for each array element 

                    origCBFs(:,:) = 0+1i*0;
                    if (true)
                        fprintf('*** debug : Number of initially generated CBFs = %d\n',numCBFsperDomain);
                    end%if
                    ind = 0;
                    for cbf_ind = cbf_start_index:cbf_end_index
                        ind = ind+1;                            
                        origCBFs(:,ind) = Ik(domain_bot_p:domain_top_p,cbf_ind);
                    end%for

                    % Plot here the SVD spectrum (only for the last, and specified iteration)
                    % and only if in local debug mode.
                    k_plot_SVD = k_iter+1;
                    plot_SVD = false;
                    if ( plot_SVD && ((k == k_iter) || (k == k_plot_SVD)) )
                        % Plot the data.
                        MBFplotSVspectrum_tmp = Const.MBFplotSVspectrum;
                        Const.MBFplotSVspectrum = true;

                        % Make sure the data is store to an *.FCD file for further processing in POSTFEKO.
                        store_to_fcd_file_tmp = Const.store_to_fcd_file;
                        Const.store_to_fcd_file = true;
                    end%if

                    % The singular values are stored to an FCD string - 
                    fcdString = sprintf('ifb_svd_alg%d_arryEl%d',Const.IFBalg,p);
                    redCBFs = reduceMBFset(Const,origCBFs,fcdString);

                    % Reset the SVD spectrum plotting and FCD file storage if enabled above.
                    if (plot_SVD && ((k == k_iter) || (k == k_plot_SVD)) )
                        Const.MBFplotSVspectrum = MBFplotSVspectrum_tmp;
                        Const.store_to_fcd_file = store_to_fcd_file_tmp;
                    end%if

                    % Increment the total number of CBFs (all domains), so that we can allocate space for Zred and
                    % Vred correctly below.
                    totRedCBFs = totRedCBFs + size(redCBFs,2);
                    
                    % Retain now the orthonormalised CBFs (also update the number of orth. CBFs / element).
                    Ik_ortho(domain_bot_p:domain_top_p,1:size(redCBFs,2)) = redCBFs;
                    %num_orth_CBFs_per_element(p) = size(redCBFs,2);
                    
                    if (true)
                        fprintf('*** debug : Number of retained orthonormal CBFs = %d\n\n',size(redCBFs,2));
                    end%if                
                    %ifbmom.svdTime(solNum) = toc; % End timing
                end% for p = 1:numArrayElements

                % Calculate here space for the reduced impedance matrix
                Zred = complex(zeros(numCBFsperDomain*numArrayElements, numCBFsperDomain*numArrayElements));
                Vred = complex(zeros(numCBFsperDomain*numArrayElements, 1));

                % Create first a column augmented vector of the primary
                % CBFs using the surface-current calculated at the current iteration(s).
                IcbfsP = complex(zeros(Nloc,numCBFsperDomain));
                for p = 1:numArrayElements
                    % Get the correct basis function range for element p:
                    domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                    domain_top_p = Const.arrayElBasisFunctRange(p,2);
                    ind = 0;
                    for cbf_ind = cbf_start_index:cbf_end_index
                        ind = ind+1;
                        % -- Use orthonormalised CBFs
                        IcbfsP(:,ind) = Ik_ortho(domain_bot_p:domain_top_p,ind);                            
                        % if (Const.useMBFreduction && ~(K_back==1))
                        %     % -- Use post-SVD reduced CBFs
                        %     IcbfsP(:,ind) = Ik_ortho(domain_bot_p:domain_top_p,ind);
                        % else
                        %     % -- Use standard CBFs (not orthonormalised)
                        %     IcbfsP(:,ind) = Ik(domain_bot_p:domain_top_p,cbf_ind);
                        % end%if
                    end%for

                    % Get now the CBFs for domain q:
                    IcbfsQ = complex(zeros(Nloc,numCBFsperDomain));
                    for q = 1:numArrayElements
                         % Get the correct basis function range for element q:
                        domain_bot_q = Const.arrayElBasisFunctRange(q,1);
                        domain_top_q = Const.arrayElBasisFunctRange(q,2);
                        ind = 0;
                        for cbf_ind = cbf_start_index:cbf_end_index
                            ind = ind+1;
                            % -- Use orthonormalised CBFs
                            IcbfsQ(:,ind) = Ik_ortho(domain_bot_q:domain_top_q,ind);
                            % if (Const.useMBFreduction&& ~(K_back==1))
                            %     % -- Use post-SVD reduced CBFs
                            %     IcbfsQ(:,ind) = Ik_ortho(domain_bot_q:domain_top_q,ind);
                            % else
                            %     % -- Use standard CBFs (not orthonormalised)
                            %     IcbfsQ(:,ind) = Ik(domain_bot_q:domain_top_q,cbf_ind);
                            % end%if
                        end%for

                        % --------------------------------------------------------------------------
                        % Calculate the correct offset for storing the entries in the reduced matrix
                        % equation (Zpq)

                        % === Range of P (rows, i.e. testing functions)
                        numMBFsP = 0;
                        for domain = 1:(p-1)
                            numMBFsP = numMBFsP + numCBFsperDomain;
                        end%for
                        PindxStart = numMBFsP+1;
                        PindxEnd   = (PindxStart - 1) + numCBFsperDomain;
                        
                        % === Range of Q (columns, i.e. basis functions)
                        numMBFsQ = 0;
                        for domain = 1:(q-1)
                            numMBFsQ = numMBFsQ + numCBFsperDomain;
                        end%for
                        QindxStart = numMBFsQ+1;
                        QindxEnd   = (QindxStart - 1) + numCBFsperDomain;

                        % -- Use newly calculated CBFs for the testing vectors
                        % Calculate the coupling matrices and Vcoupl
                        [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                        Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);
                        Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                        Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;
                    
                        % % -- Use newly calculated CBFs for the testing vectors
                        % %    Calculate the coupling matrices, Zcoupl (or Ucoupl and Vcoupl if the 
                        % %    ACA is used)
                        % [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                        % Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);
                        
                        % Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                        % Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;

                    end%for q = 1:numArrayElements
                    
                end%for p = 1, numArrayElements

                % LU-decomposition of the Z-matrix (now included here in this time)
                % to solve now for the unknown Beta coefficients, using Zred and Vred
                [Lred,Ured] = lu(Zred);
                b = Lred\Vred;
                beta_k = Ured\b;  % A one-dimensional array of length : [numCBFsperDomain * numArrayElements x 1]

                if (false)
                    % Some debug output here.
                    fprintf('*** debug : size(Zred) =(%d,%d)\n',size(Zred,1),size(Zred,2));
                    fprintf('*** debug : size(Vred) =(%d,%d)\n',size(Vred,1),size(Vred,2));
                    Zred
                    Vred
                    beta_k
                    %stop
                end%if

                % Now that we have the beta_k coefficients, we can apply them to the solution to calculate the surface
                % current at the present iteration.

                for n=1:numArrayElements % To calculate the current on domain p - as noted in [8] (used "n" here just 
                                         % for consistency witht the above code)
                    % Get the correct basis function range for element n:
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    % Extract Znn^(-1) = that of array element (with all elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;

                    for m=1:numArrayElements

                        % Add the pth coupling term for array element n using the (p-1)th coupling
                        % terms from the other array elements (discarding offcourse m=n)
                        domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                        domain_top_m = Const.arrayElBasisFunctRange(m,2);

                        % Extract first Znm
                        [Znm, Unm, Vnm] = ...
                            calcZmn(Const,zMatrices,n,m,[domain_bot_n:domain_top_n], ...
                                                        [domain_bot_m:domain_top_m]);

                        if (m~=n)
                            % Calculate the correct beta_k_offset for the array element, i.e. where its beta coefficients
                            % are stored in the beta_k vector (remember, for this algorithm we have the same number of
                            % CBFs per element - unlike algorithm 13 where it can vary - see also comment there
                            beta_k_offset = 0;
                            for array_el = 1:(m-1)
                                %beta_k_offset = beta_k_offset + num_orth_CBFs_per_element(array_el);
                                beta_k_offset = beta_k_offset + numCBFsperDomain;
                            end%for

                            ind = 0;
                            for cbf_ind = cbf_start_index:cbf_end_index
                                ind = ind+1;
                                % We use the orthonormalised CBFs
                                Ik_tmp = Ik_ortho(domain_bot_m:domain_top_m,ind);
                                % Calculate now the correct beta index (see also comment above)
                                beta_k_index = beta_k_offset + cbf_ind;

                                if (~Const.useACA)
                                    if (false)
                                        % Some debug output here.
                                        fprintf('*** debug : size(Znm) = [%d,%d]\n',size(Znm,1),size(Znm,2));
                                        fprintf('*** debug : size(Ik_tmp) = [%d,%d]\n',size(Ik_tmp,1),size(Ik_tmp,2));
                                        fprintf('*** debug : size(beta_k) = [%d,%d]\n',size(beta_k,1),size(beta_k,2));
                                        fprintf('*** debug : beta_k_offset = %d\n',beta_k_offset);
                                        fprintf('*** debug : beta_k_index  = %d\n',beta_k_index);
                                    end%if
                                    b = Ldom_n\(beta_k(beta_k_index,1)*Znm*Ik_tmp);
                                else
                                    % By using the ACA we can accelerate the matrix-vector product
                                    % By replacing Znm ~ Unm*Vnm
                                    b = Ldom_n\(beta_k(beta_k_index,1)*Unm*Vnm*Ik_tmp);
                                end%if

                                % Update surface-current at the current iteration.
                                Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                            end%for cbf_ind = cbf_start_index:cbf_end_index

                        end%if (m~=n)
                    end%for m=1:numArrayElements

                    Vdom_n = yVectors.values(domain_bot_n:domain_top_n,1);
                    b = Ldom_n\Vdom_n;
                    Idom_n = Udom_n\b;

                    % Add now this pth contribution to the current on element n
                    ifbmom.Isol(domain_bot_n:domain_top_n,1) = Idom_n - Ik(domain_bot_n:domain_top_n,k);
                    Ik(domain_bot_n:domain_top_n,k) = ifbmom.Isol(domain_bot_n:domain_top_n,1);
                end%for n=1:numArrayElements

                % Stop pre-computation timing
                ifbmom.solTime = ifbmom.solTime + toc;

                if (Const.IFB_debug>=1)
                    ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
                    ifbmom.iterTiming(k) = ifbmom.solTime;
                    if (local_debug_flag)
                        % Some debug output here.
                        fprintf('*** debug : ifbmom.relIterError(k) = %f %%\n',ifbmom.relIterError(k));
                    end%if
                end%if

                % Check whether the rel. error is <= eps_percent - if so, then exit the loop here
                if ((Const.IFB_iterations == -1) && (ifbmom.relIterError(k) <= eps_percent))
                    k_converged = k;
                    break;
                end%if

            end%for k=2:k_iter

        end%if

    % =================================================================================================================
    % 2016-12-17: This is the same as Alg. 14, but with the Orthonormalisation done per element (and not globally) and
    % by using the Gram-Schmidt orthonormalisation (and not the SVD)
    elseif (Const.IFBalg == 15)
    
        if (Const.useMBFreduction)
            message(Const, 'IFB Algorithm 15 cannot be used with SVD reduction');
            error('IFB Algorithm 14 cannot be used with SVD reduction');
        end%if

        if (Const.IFB_CBFs==1)
            message(Const, 'IFB Algorithm 14 can only be used when all previous iterations are used as CBFs');
            error('IFB Algorithm 14 can only be used when all previous iterations are used as CBFs');
        end%if

        if (Const.IFB_CBFs==0)
            message(Const, 'IFB Algorithm 14 cannot be used for a Jacobi only simulation');
            error('IFB Algorithm 14 cannot be used for a Jacobi only simulation');
        end%if

        if (calculateRelativeResiduum)
            message(Const, 'NOTE: Calculating relative residuum (execution times might be slow)');
        end%if        

        % Start timing for this algorithm
        tic

        % Store all the current values at the various iterations.
        Ik = complex(zeros(Ntot,k_iter));

        % Make sure we have enough space to store the orthonormalised eigencurrents for each array element.
        %Ik_ortho = complex(zeros(Ntot,K_back));
        Ik_ortho = complex(zeros(Ntot,k_iter));
        
        % Initialise the solution vector
        ifbmom.Isol(:,1) = 0.0d0;

        % Loop now over all the array elements and calculate the primary MBFs (starting vectors, k=1) for each
        for n=1:numArrayElements
            domain_bot_n = Const.arrayElBasisFunctRange(n,1);
            domain_top_n = Const.arrayElBasisFunctRange(n,2);
            VdomB = yVectors.values(domain_bot_n:domain_top_n,1);
            b = LdomB\VdomB;
            IdomB_0 = UdomB\b;
            Ik(domain_bot_n:domain_top_n,1) = IdomB_0;  % Domain B (array element) - primary MBF
            Ik_ortho(domain_bot_n:domain_top_n,1) = IdomB_0/sqrt((IdomB_0)'*IdomB_0); % GS ortho.
        end%for

        ifbmom.Isol(:,1) = +Ik(:,1);   % Start with k=1. Strictly speaking we should start at 0,
                                       % but MATLAB does not allow this
        % Stop pre-computation timing
        ifbmom.solTime = ifbmom.solTime + toc;

        % Note: Zero the timing here, as all the methods start off with the same solution (equal to the
        % isolated case).
        ifbmom.solTime = 0; % Zero the value here.
                                       
        if (Const.IFB_debug>=1)
            ifbmom.relIterError = zeros(1,k_iter);
            ifbmom.iterTiming = zeros(1,k_iter);

            % Also profile now the reduced matrix setup and the orthonormalisation step
            ifbmom.orthonormTiming = zeros(1,k_iter);
            ifbmom.cbfmTiming = zeros(1,k_iter);

            % Added now the relative residuum norm here also (only if enabled, as this is costly)
            if (calculateRelativeResiduum)
                ifbmom.relResError = zeros(1,k_iter);
            end

            ifbmom.relIterError(1) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
            if (calculateRelativeResiduum)
                ifbmom.relResError(1) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
            end%if

            % Initialise the timing results
            ifbmom.iterTiming(1) = ifbmom.solTime;
            ifbmom.orthonormTiming(1) = 0;
            ifbmom.cbfmTiming(1) = 0;
        end%if

        % Add secondary MBFs for each of the array elements if the number of iterations are greater
        % than 1.
        if (k_iter > 1)
            for k=2:k_iter % Iterations k= 2, ...
                
                % Time now each iteration. NOTE: We do not time the iteration as a whole now, but rather the individual
                % orthonormalization and CBFM solution steps in order to obtain where the expensive part lies (see below).
                % tic
                
                message(Const,sprintf('\nIteration k = %d',k));

                % =====================================================================================================
                % Calculate the Beta_k values here. - See [8], Alg. 3 - other than Alg. 1 and 2, the Beta coefficients
                % are calculated beforehand here.
                % =====================================================================================================

                if (false)
                    % -- Version 1 : Use 1 2 3 .... N
                    %                        {------} these CBFs
                    %                         K_back

                    % For this algorithm, we use only the latest current at the present iteration as CBF
                    if (k <= K_back)
                        cbf_start_index = 1;
                    else
                        cbf_start_index = k - K_back; % + 1; --> not correct to (+1).
                    end
                    cbf_end_index = k-1;

                else
                    % -- Version 1 : Use 1 2 3 .... N
                    %                   {------} these CBFs ... Better, as the initial ones are more orthogonal
                    %                     K_back

                    % For this algorithm, we always start from CBF index 1
                    cbf_start_index = 1;
                    if (k <= K_back)
                        cbf_end_index = k-1;
                    else
                        cbf_end_index = K_back; % + 1; --> not correct to (+1).
                    end
                end%if

                % Update now the number of CBFs that we will use.
                numCBFsperDomain = cbf_end_index - cbf_start_index + 1;
                if (true)
                    fprintf('*** debug : Number of CBFs per domain = %d\n',numCBFsperDomain);
                end%if

                % ----------------------------------------------------------------------------------------------------
                % Apply now an orthonormalisation scheme to improve the quality of the CBFs by using the Gram-Schmidt
                % ----------------------------------------------------------------------------------------------------

                % Start timing the orthonormalisation
                tic

                % Note, we start from iteration 3 and onwards here
                if (k > 2)

                    if (true)
                        fprintf('*** debug : Orthonormalising CBFs using Gram-Schmidt\n');
                    end%if

                    % Loop over all the elements and reduce the CBF set of each individually
                    for p = 1:numArrayElements
                        % Get the correct basis function range for element p:
                        domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                        domain_top_p = Const.arrayElBasisFunctRange(p,2);

                        % Note: we only update the orthogonal list with the latest CBF
                        Ik_ortho(domain_bot_p:domain_top_p,k-1) = Ik(domain_bot_p:domain_top_p,k-1);

                        for j=1:k-2 % TO-DO: Check here ....
                            Ik_ortho(domain_bot_p:domain_top_p,k-1) = Ik_ortho(domain_bot_p:domain_top_p,k-1) - ...
                            ( Ik_ortho(domain_bot_p:domain_top_p,k-1)'*Ik_ortho(domain_bot_p:domain_top_p,j) )*Ik_ortho(domain_bot_p:domain_top_p,j);
                        end%for j=1:k-2

                        % Normalise the new CBF
                        Ik_ortho(domain_bot_p:domain_top_p,k-1) = Ik_ortho(domain_bot_p:domain_top_p,k-1)/sqrt(Ik_ortho(domain_bot_p:domain_top_p,k-1)'*...
                            Ik_ortho(domain_bot_p:domain_top_p,k-1));
                    end% for p = 1:numArrayElements
                end % if (k > 2)

                % Stop the orthonormalisation timing
                ifbmom.orthonormTiming(k) = ifbmom.orthonormTiming(k-1) + toc;

                % Start timing the CBFM solution
                tic

                % Calculate here space for the reduced impedance matrix
                Zred = complex(zeros(numCBFsperDomain*numArrayElements, numCBFsperDomain*numArrayElements));
                Vred = complex(zeros(numCBFsperDomain*numArrayElements, 1));

                % Create first a column augmented column matrix of the primary CBFs
                % using the surface-current calculated at the current iteration(s).
                IcbfsP = complex(zeros(Nloc,numCBFsperDomain));
                for p = 1:numArrayElements
                    % Get the correct basis function range for element p:
                    domain_bot_p = Const.arrayElBasisFunctRange(p,1);
                    domain_top_p = Const.arrayElBasisFunctRange(p,2);
                    ind = 0;
                    for cbf_ind = cbf_start_index:cbf_end_index
                        ind = ind+1;
                        % -- Use orthonormalised CBFs
                        IcbfsP(:,ind) = Ik_ortho(domain_bot_p:domain_top_p,ind);                            
                        % if (Const.useMBFreduction && ~(K_back==1))
                        %     % -- Use post-SVD reduced CBFs
                        %     IcbfsP(:,ind) = Ik_ortho(domain_bot_p:domain_top_p,ind);
                        % else
                        %     % -- Use standard CBFs (not orthonormalised)
                        %     IcbfsP(:,ind) = Ik(domain_bot_p:domain_top_p,cbf_ind);
                        % end%if
                    end%for

                    % Get now the CBFs for domain q:
                    IcbfsQ = complex(zeros(Nloc,numCBFsperDomain));
                    for q = 1:numArrayElements
                         % Get the correct basis function range for element q:
                        domain_bot_q = Const.arrayElBasisFunctRange(q,1);
                        domain_top_q = Const.arrayElBasisFunctRange(q,2);
                        ind = 0;
                        for cbf_ind = cbf_start_index:cbf_end_index
                            ind = ind+1;
                            % -- Use orthonormalised CBFs
                            IcbfsQ(:,ind) = Ik_ortho(domain_bot_q:domain_top_q,ind);
                            % if (Const.useMBFreduction&& ~(K_back==1))
                            %     % -- Use post-SVD reduced CBFs
                            %     IcbfsQ(:,ind) = Ik_ortho(domain_bot_q:domain_top_q,ind);
                            % else
                            %     % -- Use standard CBFs (not orthonormalised)
                            %     IcbfsQ(:,ind) = Ik(domain_bot_q:domain_top_q,cbf_ind);
                            % end%if
                        end%for

                        % --------------------------------------------------------------------------
                        % Calculate the correct offset for storing the entries in the reduced matrix
                        % equation (Zpq)

                        % === Range of P (rows, i.e. testing functions)
                        numMBFsP = 0;
                        for domain = 1:(p-1)
                            numMBFsP = numMBFsP + numCBFsperDomain;
                        end%for
                        PindxStart = numMBFsP+1;
                        PindxEnd   = (PindxStart - 1) + numCBFsperDomain;
                        
                        % === Range of Q (columns, i.e. basis functions)
                        numMBFsQ = 0;
                        for domain = 1:(q-1)
                            numMBFsQ = numMBFsQ + numCBFsperDomain;
                        end%for
                        QindxStart = numMBFsQ+1;
                        QindxEnd   = (QindxStart - 1) + numCBFsperDomain;

                        % See issue FEKDDM-3.2: Add now ACA support for the coupling matrix assebly here (only for off-diagonal terms)
                        if (p~=q)                            
                            [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);
                        else
                            % DO not use ACA for diagonal terms
                            useACAtmp = Const.useACA;
                            Const.useACA = 0;
                            [Zcoupl, Ucoupl, Vcoupl] = calcZmn(Const,zMatrices,p,q,[domain_bot_p:domain_top_p],[domain_bot_q:domain_top_q]);                            
                            %Zcoupl = zMatrices.values((p-1)*Ndom+1:p*Ndom,(q-1)*Ndom+1:q*Ndom);
                            Const.useACA = useACAtmp;
                        end
                        
                        % Calculate Vcoupl
                        Ycoupl = yVectors.values([domain_bot_p:domain_top_p],1);

                        % ------------------------------------------------------------------------------
                        % Block assignment for submatrix in Zred.
                        % See issue FEKDDM-3.2 on Basecamp (and also above comment): The ACA can also be used for Zcoupl
                        if ((Const.ACAalg < 3) || (p==q))                            
                            Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Zcoupl*IcbfsQ;
                        else
                            Zred(PindxStart:PindxEnd,QindxStart:QindxEnd) = ((IcbfsP).')*Ucoupl*Vcoupl*IcbfsQ;
                        end%if


                        Vred(PindxStart:PindxEnd) = ((IcbfsP).')*Ycoupl;
                    end%for q = 1:numArrayElements
                    
                end%for p = 1, numArrayElements

                % calculate the memory usage here of Zred:
                ifbmom.memUsageZred = byteSize(Zred);

                % LU-decomposition of the Z-matrix (now included here in this time)
                % to solve now for the unknown Beta coefficients, using Zred and Vred
                [Lred,Ured] = lu(Zred);
                b = Lred\Vred;
                beta_k = Ured\b;  % A one-dimensional array of length : [numCBFsperDomain * numArrayElements x 1]

                if (false)
                    % Some debug output here.
                    fprintf('*** debug : size(Zred) =(%d,%d)\n',size(Zred,1),size(Zred,2));
                    fprintf('*** debug : size(Vred) =(%d,%d)\n',size(Vred,1),size(Vred,2));
                    Zred
                    Vred
                    beta_k
                    %stop
                end%if

                % Now that we have the beta_k coefficients, we can apply them to the solution to calculate the surface
                % current at the present iteration.

                % Switch off ACA calculation here for the mat-vec below
                useACAtmp = Const.useACA;
                Const.useACA = 0;

                for n=1:numArrayElements % To calculate the current on domain p - as noted in [8] (used "n" here just 
                                         % for consistency witht the above code)
                    % Get the correct basis function range for element n:
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    % Extract Znn^(-1) = that of array element (with all elements assumed identical).
                    Udom_n = UdomB;
                    Ldom_n = LdomB;

                    for m=1:numArrayElements

                        % Add the pth coupling term for array element n using the (p-1)th coupling
                        % terms from the other array elements (discarding offcourse m=n)
                        domain_bot_m = Const.arrayElBasisFunctRange(m,1);
                        domain_top_m = Const.arrayElBasisFunctRange(m,2);                

                        % Extract first Znm
                        [Znm, Unm, Vnm] = ...
                            calcZmn(Const,zMatrices,n,m,[domain_bot_n:domain_top_n], ...
                                                        [domain_bot_m:domain_top_m]);                        

                        if (m~=n)
                            % Calculate the correct beta_k_offset for the array element, i.e. where its beta coefficients
                            % are stored in the beta_k vector (remember, for this algorithm we have the same number of
                            % CBFs per element - unlike algorithm 13 where it can vary - see also comment there
                            beta_k_offset = 0;
                            for array_el = 1:(m-1)
                                %beta_k_offset = beta_k_offset + num_orth_CBFs_per_element(array_el);
                                beta_k_offset = beta_k_offset + numCBFsperDomain;
                            end%for

                            ind = 0;
                            for cbf_ind = cbf_start_index:cbf_end_index
                                ind = ind+1;
                                % We use the orthonormalised CBFs
                                Ik_tmp = Ik_ortho(domain_bot_m:domain_top_m,ind);
                                % Calculate now the correct beta index (see also comment above)
                                beta_k_index = beta_k_offset + cbf_ind;
                                
                                if (~Const.useACA)
                                    if (false)
                                        % Some debug output here.
                                        fprintf('*** debug : size(Znm) = [%d,%d]\n',size(Znm,1),size(Znm,2));
                                        fprintf('*** debug : size(Ik_tmp) = [%d,%d]\n',size(Ik_tmp,1),size(Ik_tmp,2));
                                        fprintf('*** debug : size(beta_k) = [%d,%d]\n',size(beta_k,1),size(beta_k,2));
                                        fprintf('*** debug : beta_k_offset = %d\n',beta_k_offset);
                                        fprintf('*** debug : beta_k_index  = %d\n',beta_k_index);
                                    end%if
                                    b = Ldom_n\(beta_k(beta_k_index,1)*Znm*Ik_tmp);
                                else
                                    % By using the ACA we can accelerate the matrix-vector product
                                    % By replacing Znm ~ Unm*Vnm
                                    b = Ldom_n\(beta_k(beta_k_index,1)*Unm*Vnm*Ik_tmp);
                                end%if

                                % Update surface-current at the current iteration.
                                Ik(domain_bot_n:domain_top_n,k) = Ik(domain_bot_n:domain_top_n,k) + Udom_n\b;
                            end%for cbf_ind = cbf_start_index:cbf_end_index

                        end%if (m~=n)
                    end%for m=1:numArrayElements

                    Vdom_n = yVectors.values(domain_bot_n:domain_top_n,1);
                    b = Ldom_n\Vdom_n;
                    Idom_n = Udom_n\b;

                    % Add now this pth contribution to the current on element n
                    ifbmom.Isol(domain_bot_n:domain_top_n,1) = Idom_n - Ik(domain_bot_n:domain_top_n,k);
                    Ik(domain_bot_n:domain_top_n,k) = ifbmom.Isol(domain_bot_n:domain_top_n,1);
                end%for n=1:numArrayElements

                % See above comment - switch the ACA back to its original status
                Const.useACA = useACAtmp;  

                % Stop the CBFM timing
                ifbmom.cbfmTiming(k) = ifbmom.cbfmTiming(k-1) + toc;

                if (Const.IFB_debug>=1)
                    ifbmom.relIterError(k) = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);
                    % Calculate the relative residual error here (only if enabled, as this is expensive)
                    if (calculateRelativeResiduum)
                        ifbmom.relResError(k) = (calculateErrorNormPercentage(yVectors.values, zMatrices.values*ifbmom.Isol))/100.0;
                    end%if
                    % Calculate the total time - sum of the Orthonormalisation and CBFM timing
                    ifbmom.iterTiming(k) = ifbmom.orthonormTiming(k) + ifbmom.cbfmTiming(k);
                    if (local_debug_flag)
                        % Some debug output here.
                        fprintf('*** debug : ifbmom.relIterError(k) = %f %%\n',ifbmom.relIterError(k));
                    end%if
                end%if

                % Check whether the rel. error is <= eps_percent - if so, then exit the loop here
                if ((Const.IFB_iterations == -1) && (ifbmom.relIterError(k) <= eps_percent))
                    k_converged = k;
                    break;
                end%if

            end%for k=2:k_iter

        end%if

        % Make sure we record the total time again:
        ifbmom.solTime = ifbmom.iterTiming(k_iter);

    end%if


    % ----------------------------------------------------------------------------------------------
    % Restore the following configuration parameters if Const.IFB_alg=6 or 8 (i.e. when the DGFM
    % is used)
    if ((Const.IFBalg == 6) || (Const.IFBalg == 8) )
        Const.runDGFMfromIFBMoMsolver = false;
        Const.DGFMweightVectorCalcScheme = DGFMweightVectorCalcScheme_tmp;
        Const.useDGFMmethod = useDGFMmethod_tmp;
        Const.runJACKITfromDGFM = runJACKITfromDGFM_tmp;
        % Also restore the excitation vector
        yVectors = yVectors_tmp;
    end%if

    % ==============================================================================================
    % TO-DO: This needs more work!
    % Calculate the memory usages
    if (Const.IFBalg == 7)
        % -- Jacobi method
        %    size(Zpp0^(-1)) + size(Zpq)
        ifbmom.memUsage = byteSize(complex(zeros(size(ZdomB,1)*2,size(ZdomB,2)*1)));
    elseif(Const.IFBalg == 9)
        % -- CBFM-enhnaced Jacobi Method (Alg. 1)
        ifbmom.memUsage = byteSize(complex(zeros(size(ZdomB,1)*2,size(ZdomB,2)*1)));
    
    % NO LONGER IN USE BELOW -- NOT CORRECTLY ACCOUNTED FOR
    % elseif ((Const.IFBalg == 11)||(Const.IFBalg == 12)||(Const.IFBalg == 13)||(Const.IFBalg == 14))
    %     % -- CBFM-enhnaced Jacobi Method (Alg. 11 - 14)
        
    %     ifbmom.memUsage = byteSize(complex(zeros(size(ZdomB,1)*2,size(ZdomB,2)*1)));
    %     %ifbmom.memUsage = byteSize(Z0); 
    %     %ifbmom.memUsage = byteSize(V0);
    
    elseif ((Const.IFBalg == 13) || (Const.IFBalg == 15))
        % -- CBFM-enhnaced Jacobi Method (Alg. 15) : Need to add also the CBFM red matrix
        
        ifbmom.memUsage = byteSize(complex(zeros(size(ZdomB,1)*2,size(ZdomB,2)*1)));
        %ifbmom.memUsage = byteSize(Z0); 
        %ifbmom.memUsage = byteSize(V0);
    else
        % -- Just set it equal to the MoM
        ifbmom.memUsage = byteSize(zMatrices.values);% + byteSize(Zsm);
    end%if
    % ==============================================================================================

    % Compare the solution obtained here, with that obtained by FEKO that was stored in xVectors.values
    if ( (Const.IFBalg~=5) && (Const.IFBalg~=6) && (Const.IFBalg~=7) && (Const.IFBalg~=8) && (Const.IFBalg~=9) ...
        && (Const.IFBalg~=10) && (Const.IFBalg~=11) && (Const.IFBalg~=12) && (Const.IFBalg~=13) && (Const.IFBalg~=14) ...
        && (Const.IFBalg~=15) )
        ifbmom.relErrorDomA = calculateErrorNormPercentage(xVectors.values(1:NdomA), ifbmom.Isol(1:NdomA));
    end%if
    ifbmom.relError = calculateErrorNormPercentage(xVectors.values, ifbmom.Isol);

    % TO-DO: Danie, support the following:
    % Write the NGF-en. DGFM solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested
    % if (isequal(and([1 0 0 0],Const.FEKDDMwriteFEKOstrfile),[1 0 0 0]))
    %     writeSolToFile(Const, ngfdgfm);
    % end%if

    % Write the IFB-MoM solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested
    if (isequal(and([0 1 0 0 0 0 0],Const.FEKDDMwriteFEKOstrfile),[0 1 0 0 0 0 0]))
        writeSolToFile(Const, ifbmom);
    end%if

    message(Const,sprintf('Finished IFB-MoM (total) solver in %f sec.',ifbmom.solTime));
    %message(Const,sprintf('(Setup time %f sec.)',ifbmom.solTime));

    message(Const,sprintf('Memory usage of IFB-MoM %s',ifbmom.memUsage));

    % 2017-05-16: Also print out here the memory usage for storing Zred (CBFM red. matrix) for algorithms 13 and 15
    if ((Const.IFBalg == 13) || (Const.IFBalg == 15))
        message(Const,sprintf('Memory usage of IFB-MoM storing Zred %s',ifbmom.memUsageZred));
    end%if

    message(Const,sprintf('Rel. error norm. compared to FEKO sol. %f percent:',ifbmom.relError));
    if ((Const.IFBalg~=5) && (Const.IFBalg~=6) && (Const.IFBalg~=7) && (Const.IFBalg~=8)&&(Const.IFBalg~=9)...
        && (Const.IFBalg~=10) && (Const.IFBalg~=11) && (Const.IFBalg~=12) && (Const.IFBalg~=13) && (Const.IFBalg~=14) ...
        && (Const.IFBalg~=15))
        message(Const,sprintf('Rel. error norm. compared to FEKO sol. - (dom. A) %f precent:',ifbmom.relErrorDomA));
    end%if

    if ((Const.IFB_iterations == -1) && ((Const.IFBalg==7) || (Const.IFBalg==9) || (Const.IFBalg==11) || ...
        (Const.IFBalg==12) || (Const.IFBalg==13) || (Const.IFBalg==14) || (Const.IFBalg==15)))
        message(Const,sprintf('Convergence after %d iterations for eps = %f percent', k_converged,eps_percent));
    end%if

    % Some additional post-processing (if IFB-MoM solver >= 1)
    if (Const.IFB_debug>=1)
        % Account here for the fact that we might have stopped after a certain number of iterations, if the 
        % threshold was reached.
        y1=ifbmom.relIterError;
        if (Const.IFB_iterations == -1)
            x1=[1:k_converged];
            y1=y1(1:k_converged);
        else
            x1=[1:k_iter];
        end%if
        
        xlab = 'Iteration (K)';
        ylab = 'Relative error norm percentage (\epsilon) [%]';
        % Build the title string here.
        if ((Const.IFBalg == 12) || (Const.IFBalg == 13) || (Const.IFBalg == 14) || (Const.IFBalg == 15))
            titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                Const.IFBalg,K_back, Const.useMBFreduction);
        else
            titleString = sprintf('IFB algorithm : %d',Const.IFBalg);
        end%if
        imgString = strcat('eps_',img_fcd_filename);
        fcdString = imgString; % Name the *.fcd name, same as the image file name.
        Const.plotSemiLogY=true;  
        plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
        Const.plotSemiLogY=false;

        if (calculateRelativeResiduum)
            % Account here for the fact that we might have stopped after a certain number of iterations, if the 
            % threshold was reached.
            y1=ifbmom.relResError;
            if (Const.IFB_iterations == -1)
                x1=[1:k_converged];
                y1=y1(1:k_converged);
            else
                x1=[1:k_iter];
            end%if
            
            xlab = 'Iteration (K)';
            ylab = 'Relative residuum (\epsilon)';
            % Build the title string here.
            if ((Const.IFBalg == 12) || (Const.IFBalg == 13) || (Const.IFBalg == 14) || (Const.IFBalg == 15))
                titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                    Const.IFBalg,K_back, Const.useMBFreduction);
            else
                titleString = sprintf('IFB algorithm : %d',Const.IFBalg);
            end%if
            imgString = strcat('relres_',img_fcd_filename);
            fcdString = imgString; % Name the *.fcd name, same as the image file name.
            Const.plotSemiLogY=true;  
            plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
            Const.plotSemiLogY=false;
        end%if

        % Plot the timing - if we are working with Alg. 12 - 15
        % 2017-05-16: Do not always write out the files here (Set a flag at the end to control this) - NOT SURE WHY THIS WAS ADDED?
        if (((Const.IFBalg == 12) || (Const.IFBalg == 13) || (Const.IFBalg == 14) || (Const.IFBalg == 15)) && true)
            

            % --------------------------
            % Plot the iteration timing
            % --------------------------

            % Account here for the fact that we might have stopped after a certain number of iterations, if the 
            % threshold was reached.
            y1=ifbmom.iterTiming;
            if (Const.IFB_iterations == -1)
                x1=[1:k_converged];
                y1=y1(1:k_converged);
            else
                x1=[1:k_iter];
            end%if
            
            % Store here the iteration timing
            iterTiming = y1;

            xlab = 'Iteration (K)';
            ylab = 'Total runtime [seconds]';
            % Build the title string here.
            titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                    Const.IFBalg,K_back, Const.useMBFreduction);
            imgString = strcat('time_',img_fcd_filename);
            fcdString = imgString; % Name the *.fcd name, same as the image file name.
            Const.plotSemiLogY=true;  
            plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
            Const.plotSemiLogY=false;

            % --------------------------
            % Plot the orthonormalisation timing (not for Jacobi method, i.e. Algorithm 12 with number of CBFs = 0)
            % --------------------------

            if (~use_Jacobi) % Only plot the orthonormalisation if the Jacobi method has not been used
                % Account here for the fact that we might have stopped after a certain number of iterations, if the 
                % threshold was reached.
                y1=ifbmom.orthonormTiming;
                if (Const.IFB_iterations == -1)
                    x1=[1:k_converged];
                    y1=y1(1:k_converged);
                else
                    x1=[1:k_iter];
                end%if

                xlab = 'Iteration (K)';
                ylab = 'Orthonormalisation runtime [seconds]';
                % Build the title string here.
                titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                        Const.IFBalg,K_back, Const.useMBFreduction);
                imgString = strcat('orthotime_',img_fcd_filename);
                fcdString = imgString; % Name the *.fcd name, same as the image file name.
                Const.plotSemiLogY=true;  
                plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
                Const.plotSemiLogY=false;

                % --------------------------
                % Plot the CBFM (+ rest of alg.) timing
                % --------------------------

                % Account here for the fact that we might have stopped after a certain number of iterations, if the 
                % threshold was reached.
                y1=ifbmom.cbfmTiming;
                if (Const.IFB_iterations == -1)
                    x1=[1:k_converged];
                    y1=y1(1:k_converged);
                else
                    x1=[1:k_iter];
                end%if

                xlab = 'Iteration (K)';
                ylab = 'CBFM + rest runtime [seconds]';
                % Build the title string here.
                titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                        Const.IFBalg,K_back, Const.useMBFreduction);
                imgString = strcat('cbfmtime_',img_fcd_filename);
                fcdString = imgString; % Name the *.fcd name, same as the image file name.
                Const.plotSemiLogY=true;  
                plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
                Const.plotSemiLogY=false;
            end%if (~use_Jacobi)

            % =======================
            % New addition: 2016-08-02:
            % -------------------------
            % Plot now the relative error norm percentage convergence as a function of the runtime:
            x1 = iterTiming;
            y1 = ifbmom.relIterError;
            xlab = 'Total runtime [seconds]';
            ylab = 'Relative error norm percentage (\epsilon) [%]';
            % Build the title string here.
            titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                    Const.IFBalg,K_back, Const.useMBFreduction);
            imgString = strcat('epsvstime_',img_fcd_filename);
            fcdString = imgString; % Name the *.fcd name, same as the image file name.
            Const.plotSemiLogY=true;  
            plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
            Const.plotSemiLogY=false;

            % =======================
            % New addition: 2017-01-08:
            % -------------------------
            % Plot now the relative residual as a function of the runtime:
            x1 = iterTiming;
            y1=ifbmom.relResError;
            xlab = 'Total runtime [seconds]';
            ylab = 'Relative residuum (\epsilon)';
            % Build the title string here.
            titleString = sprintf('IFB algorithm : %d, K-back : %d, useMBFreduction : %d',...
                    Const.IFBalg,K_back, Const.useMBFreduction);
            imgString = strcat('relresvstime_',img_fcd_filename);
            fcdString = imgString; % Name the *.fcd name, same as the image file name.
            Const.plotSemiLogY=true;  
            plotData1(Const,x1,y1,xlab,ylab,titleString,imgString,fcdString);
            Const.plotSemiLogY=false;

        end%if

        % Compute also the spectral radius of the iteration matrix:
        if (Const.IFB_debug>=2)

            % TO-DO: The following is not really valid for (Const.IFBalg == 9)
            if ((Const.IFBalg == 9) || (Const.IFBalg == 10) || (Const.IFBalg == 11) || (Const.IFBalg == 12) || ...
                (Const.IFBalg == 13) || (Const.IFBalg == 14) || (Const.IFBalg == 15))
                message(Const,sprintf('ERROR: IFB algorithm 9 not valid here'));
                error(['ERROR: IFB algorithm 9 not valid here']);
            end%if

            % diagonal part of Z and rest
            Zon=complex(zeros(Ntot,Ntot));
            % TO-DO: Danie, Alg 6 or 8 is not really valid here (as the iteration matrix is
            % replaced with an active impedance matrix equation)
            if ((Const.IFBalg == 5) || (Const.IFBalg == 6) || (Const.IFBalg == 7) || (Const.IFBalg == 8))
                for n=1:numArrayElements
                    domain_bot_n = Const.arrayElBasisFunctRange(n,1);
                    domain_top_n = Const.arrayElBasisFunctRange(n,2);
                    Zon(domain_bot_n:domain_top_n,domain_bot_n:domain_top_n) = ZdomB;
                end%for
            else
                Zon(1:NdomA,1:NdomA) = ZdomA;
                Zon(NdomA+1:Ntot,NdomA+1:Ntot) = ZdomB;
            end%if

            Zoff = zMatrices.values - Zon;

            % iteration matrix (TO-DO: Check whether we need the + or - sign here). A few tests
            % showed that it does not matter.
            T = - inv(Zon) * Zoff;

            % Spectral radius condition
            ifbmom.rho = max(abs(eig(T)));
            message(Const,sprintf('Spectral radius for IFB-MoM %f',ifbmom.rho));
            if ifbmom.rho >= 1
                message(Const,sprintf('WARNING: Convergence issue possible'));
            end%if
        end%if

    end%if