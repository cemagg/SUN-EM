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
    
    % local debug flag
    LOCAL_DEBUG = true;
    
    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running DGFM interpolation solver'));
    message_fc(Const,sprintf('  Interpolation scheme to be used: %d.',Const.useDGFMinterpolation));

    % Currently we only allow this for a single solution
    numSols = yVectors.numRhs;
    if (numSols ~= 1)
        message_fc(Const,sprintf('[dgfmInterpolate] Only a single solution configuration allowed'));
        error('[dgfmInterpolate]  Only a single solution configuration allowed');
    end%if
    
    % Check also for disconnected domains
    if (~Solver_setup.disconnected_domains)
        message_fc(Const,sprintf('[dgfmInterpolate] Only disconnected array geometries allowed'));
        error('[dgfmInterpolate] Only disconnected array geometries allowed');
    end%if
    
    % Initialisations
    dgfm_interpolate  = dgfm;    
    numArrayEls = Solver_setup.num_finite_array_elements;  % The number of array elements
    
    % Now that we have the interpolated MBF coefficients (that we obtain
    % using the same set of MBFs / element - i.e. that of the reference
    % element, which in our case is set to some element near the centre)
    % we can interpolate to extract now the MBFs on the remaining array
    % indices.
    
    tic % Start timing
        % Extrapolate over each of the MBFs
        for MBFindex = 1:dgfm.numRefMBFs
            % Which MBF index do we want to extrapolate?

            % Create the interpolant using the complex MBFs
            
            % NOTE: This is WRONG below ... we need to calculate the
            % interpolant given the SAMPLES ... as it will be zero on most!
            
            F = scatteredInterpolant(dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2),...
                dgfm_interpolate.MBF_weights(:,MBFindex));
            
            % Evaluate the interpolant on the full array lattice
            Vq = F(dgfm_interpolate.array_XY(:,1),dgfm_interpolate.array_XY(:,2));
            
            % Calculate now the interpolated MBF weights.
            dgfm_interpolate.MBF_weights_interpolated(:,MBFindex) = Vq;
            
        end%for MBFindex = 1:dgfm.numRefMBFs
            
    dgfm.interpolant_calculation = toc; % End timing
        
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

    % Once we have the interpolated MBF coefficients for each of the reference
    % MBFs, we can generate now the solution on the missing array elements.
    % For each of the array elements, we need to extract a set of MBFs, with their
    % respective coefficients. The set of reduced MBFs have already been setup for each 
    % of the elements in the runMBFgenerator routine.
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
    end%for