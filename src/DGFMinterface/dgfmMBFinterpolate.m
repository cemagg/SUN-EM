function [dgfm_interpolate] = dgfmMBFinterpolate(Const, Solver_setup, yVectors, dgfm, mbfs)
    %runDGFMsolver
    %   Usage:
    %       [dgfm] = runDGFMsolver(Const, Solver_setup, zMatrices, yVectors, xVectors, mbfs, ngf)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details        
    %       yVectors
    %           The Yrhs-vector data
    %       dgfm
    %           The DGFM data (e.g. active impedance matrices) associated
    %           with each of the array elements.
    %       mbfs
    %           The MBFs (prim & sec. / reduced set) generated in runMBFgenerator
    %
    %   Output Arguments:
    %       dgfm
    %           Structs containing DGFM solution calculated after interpolation and timing data
    %
    %   Description:
    %       Given as input the DGFM active impedance matrices as well as a
    %       MBF set for each array element (typicall a SVD reduced set), we
    %       first obtain MBF coefficients for all the elements. We then
    %       select only a few and see whether the others can be calculated
    %       using complex interpolation, similar to that of Ngoy and Dirk
    %       with the CBFPs - see [1].
    %
    %   References:
    %   [1] Dirk de Villers, Ngoy, CBFP TAP article (TO-DO: Add reference)

    narginchk(5,5);
    
    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running DGFM interpolation solver'));
    message_fc(Const,sprintf('  Number of DGFM basis functions: %d.',Solver_setup.num_mom_basis_functions));

    % Currently we only allow this for a single solution
    numSols = yVectors.numRhs;
    if (numSols ~= 1)
        message_fc(Const,sprintf('[dgfmInterpolate] Only a single solution configuration allowed'));
        error('[dgfmInterpolate]  Only a single solution configuration allowed');
    end%if
    
    % TO-DO: Check also for disconnected domains
    if (~Solver_setup.disconnected_domains)
        message_fc(Const,sprintf('[dgfmInterpolate] Only disconnected array geometries allowed'));
        error('[dgfmInterpolate] Only disconnected array geometries allowed');
    end%if
    
    % Initialisations
    dgfm_interpolate  = dgfm;
    
    numArrayEls = Solver_setup.num_finite_array_elements;  % The number of array elements
    
    % Allocate some space for the MBF coefficients associated with each
    % array element. TO-DO: We need to rather keep the number of MBFs /
    % domain consistent. As an initial test, we keep the MBFs for the 
    % first element in the array.
    ref_domain = 1;
    numMBFs = mbfs.numRedMBFs(ref_domain,1);
    ref_domain_basis_functions = Solver_setup.rwg_basis_functions_domains{ref_domain};
    MBFs = mbfs.RedIsol(ref_domain_basis_functions,1:numMBFs,ref_domain,1);
            
    dgfm_interpolate.MBF_weights = complex(zeros(numArrayEls,numMBFs));
    
    % For each of the arrays, we need to extract a set of MBFs, with their
    % respective coefficients. The set of reduced MBFs have already been
    % setup for each of the elements in the runMBFgenerator routine.
    for m=1:numArrayEls
        domain_basis_functions = Solver_setup.rwg_basis_functions_domains{m};
        
        % Extract the MBFs associated with this domain.
        %MBFs = mbfs.RedIsol(domain_basis_functions,1:numMBFs,m,1);
        
        % Now that we have the MBFs, we can use the active impedance
        % matrix obtained from the i-DGFM solution to calculate the complex
        % coefficients should we expand the current on the mth domain using
        % the set of reduced MBFs. Note, the index is usually calculated as
        % follows: index = solNum + (freq-1)*numRHSperFreq;
        % For now however, we assume only a single solution number and a 
        % single frequency
        Zact = dgfm.Zact(domain_basis_functions,domain_basis_functions,1);
        
        % Now we can find the MBF weights by setting up a reduced impedance
        % matrix
        Zred = (MBFs)' * Zact * MBFs;
        Vrwg = yVectors.values(domain_basis_functions,1);
        Vred = (MBFs)' * Vrwg;
        MBF_weights = Zred \ Vred;
        
        % Store these MBF coefficients now for later processing.
        dgfm_interpolate.MBF_weights(m,:) = MBF_weights;
        
        % Compare now this solution against that obtained from the DGFM.
        if (true)
            % First build again Isol using the MBFs
            Isol_m = MBFs * MBF_weights;
            relError = calculateErrorNormPercentage(dgfm.Isol(domain_basis_functions,1), Isol_m);
            message_fc(Const,sprintf('Rel. error norm. for element. %d compared to DGFM sol. %f percent',m, relError));
        end
    end%for
    
    % Now that we have the interpolated MBF coefficients (that we obtain
    % using the same set of MBFs / element - i.e. that of the reference
    % element, which in our case is set to some element near the centre).
    
    % Extract the array positions
    dgfm_interpolate.array_XY = zeros(numArrayEls,2);    
    [dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2)] = extract_array_element_positions(Const, 'array_layout.xml');
    
    figure
    plot(dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2), 'o', 'MarkerFaceColor','k','MarkerSize',10);
    
    % Plot now the real part of the first MBF coefficient on this grid
    %VqReal = real(dgfm_interpolate.MBF_weights(:,1));
    %figure
    %plot3(dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2),VqReal,'o');
    
    % Select now only a few samples of the actual array layout at which the
    % MBF coefficients are known (the other points are unknown and have to
    % be calculated)
    
    
    % Which MBF index do we want to extrapolate?
    MBFindex = 2;
    
    tic
        % Create the interpolant using the complex MBFs
        F = scatteredInterpolant(dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2),...
            dgfm_interpolate.MBF_weights(:,MBFindex));
    toc;
        
    % Just sample this function on a regular grid for visualisation purposes.
    [Xq,Yq] = meshgrid(-9:0.5:9);
    Vq = F(Xq,Yq);
    
    % Plot the real part of the first MBF coefficient
    VqReal = real(Vq);
    figure
    hold on;
    surf(Xq,Yq,VqReal);
    plot3(dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2), real(dgfm_interpolate.MBF_weights(:,MBFindex)),'o');
    xlabel('X');
    ylabel('Y');
    zlabel('Real Value - \beta');
    title(sprintf('Real Component of MBF coeff. %d', MBFindex));
   
%     % Plot the imaginary part of the first MBF coefficient
%     VqImag = imag(Vq);
%     figure
%     surf(Xq,Yq,VqImag);
%     xlabel('X');
%     ylabel('Y');
%     zlabel('Imaginary Value - \beta');
%     title(sprintf('Imaginary Component of MBF coeff. %d',MBFindex));