function Const = sunem_init(Const, yVectors)
    %fc_init
    %   Usage:
    %           fc_init(Const, yVectors)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging settings of which solver to run,
    %           as well as initial basis function setup
    %       yVectors
    %           The Yrhs-vector data
    %
    %   Output Arguments:
    %       Const
    %           Struct containing updated configuration parameters
    %
    %   Description:
    %       Initialises the FEKO-Connect environment. Useful to make sure that certain
    %       parameters have default values, so that if they are not set in the driver
    %       script, the execution can still continue. These configuration parameters
    %       are typically passed in the Const struct (setup in the driver routine).
    %   =======================
    %   Written by Danie Ludick on 2018.01.01
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    %   =======================     

    narginchk(2,2);

    message_fc(Const,' ');
    message_fc(Const,...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Initialising the configuration setup'));

    % =======================
    % Specify a subset of solution configurations in the driver routine. 
    % Check here whether this is set and if this is in the correct range, 
    % i.e. corresponds to numSols:
    numSols = yVectors.numRhs;

        % ================================    
    try
        set = false;
        if (Const.runMoMsolver)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.runMoMsolver = false;
        end
    end
    
   % ================================    
    try
        set = false;
        if (Const.runCBFMsolver)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.runCBFMsolver = false;
        end
    end

    % ================================    
    try
        set = false;
        if (Const.runJacobisolver)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.runJacobisolver = false;
        end
    end

    % ================================    
    try
        set = false;
        if (Const.runIFBMoMsolver)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.runIFBMoMsolver = false;
        end
    end

    % ================================    
    try
        set = false;
        if (Const.runDGFMsolver)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.runDGFMsolver = false;
        end
    end
    
    % ================================    
    % First check whether the variables exist:
    % -- solStart
    try
        set = false;
        if (Const.solStart)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            Const.solStart = 1;
        end
    end
    
    % -- solEnd
    try
        set = false;
        if (Const.solEnd)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            Const.solEnd = numSols;
        end
    end

    % Check that the solution settings (if set) is indeed in the correct
    % range
    if ( ((Const.solStart < 0) || (Const.solStart > numSols)) || ...
            ((Const.solEnd   < 0) || (Const.solEnd   > numSols)) || ...
            (Const.solStart > Const.solEnd) )
        % Invalid value for solStart and solEnd
        message_fc(Const,'[fc_init] Invalid value for solStart and solEnd');
        error ('[fc_init] Invalid value for solStart and solEnd');
    end

    % Use fast MoM block calculation by default if not set
    % i.e. a vectorised implementation.
    try
        set = false;
        if (Const.fastBuildMoMblock)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Default to improved (faster) algorithm
            Const.fastBuildMoMblock = true;
        end
    end
    
    % Use the Equivalent Dipole Method (EDM) for the Z-matrix calc.
    % by default not activated. Note, this is for the internal 
    % MoM matrix calculation.
    try
        set = false;
        if (Const.useEDM)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)            
            Const.useEDM = false;
        end
    end

    % ================================    
    try
        set = false;
        if (Const.parse_FEKO_out_file)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % No debug output for the IFB-MoM solver
            Const.parse_FEKO_out_file = false;
        end
    end

    % ================================    
    try
        set = false;
        if (Const.no_mutual_coupling_array)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Deactivate this option
           Const.no_mutual_coupling_array = false;
        end
    end

    % ================================    
    % Note, the following will be overwritten below
    try
        set = false;
        if (Const.domain_decomposition)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Deactivate this option
           Const.domain_decomposition = false;
        end
    end
    
    % ================================    
    try
        set = false;
        if (Const.useMBFreduction)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Deactivate this option
           Const.useMBFreduction = false;
        end
    end
    
    % ================================    
    try
        set = false;
        if (Const.MBFthreshold)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Deactivate this option
           Const.MBFthreshold = 1000;
        end
    end

    % ================================    
    try
        set = false;
        if (Const.JACKITcheckConvergence)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Deactivate this option
           Const.JACKITcheckConvergence = false;
        end
    end

    % Also do the same with the output files
    % ================================    
    try
        set = false;
        if (Const.SUNEMcbfmstrfilename)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.SUNEMcbfmstrfilename = '';
        end
    end

    % ================================    
    try
        set = false;
        if (Const.SUNEMifbmomstrfilename)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.SUNEMifbmomstrfilename = '';
        end
    end
    
    % ================================    
    try
        set = false;
        if (Const.SUNEMmomstrfilename)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.SUNEMmomstrfilename = '';
        end
    end

    % ================================    
    try
        set = false;
        if (Const.SUNEMdgfmstrfilename)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.SUNEMdgfmstrfilename = '';
        end
    end

    % ================================    
    try
        set = false;
        if (Const.SUNEMjackstrfilename)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % Solver not activated
           Const.SUNEMjackstrfilename = '';
        end
    end

    % =================================
    % Other initialisations:
    Const.useACA = false; % ACA interface available, but not yet activated. (Requires some refactoring).
    Const.ACAalg = 3; % When ACA is activated, we use the efficient algorithm (3).
    Const.store_to_fcd_file = false; % Store to external *.fcd file - useful when importing into POSTFEKO.
    Const.useCSCBFM = false;
    Const.isPhasedArray = false;

    % Deactivate for now
    Const.runJACKITfromDGFM = false;
    Const.runNGFenDGFMsolver = false;
    Const.storeZact = false;

    % =================================
    % General constants:
    Const.RAD2DEG = 180/pi;  % For converting between radians and degrees
    Const.DEG2RAD = pi/180;  % For converting between degrees and radians
    Const.EPS = 10e-6;       % Set the machine precision to be used

    Const.EPS_0 = 8.854e-12; % Permittivity of free space
    Const.MU_0=4*pi*1e-7;    % Permeability of free space
    Const.ETA_0 = sqrt(Const.MU_0/Const.EPS_0); % Wave impedance (~ 120*pi)
    Const.C0  = 1/sqrt(Const.EPS_0*Const.MU_0); % Speed of light in free space (~ 3E8 m/s)

    % For numerical integration over a triangular domain (based on DBD2011)
    Const.QUAD_PTS = 6;     % Quadrature rule (6 point is a good default)
    Const.SING = false;     % If set to true, then singularities are treated in the DBD2011 routine for
                            % filling the Z matrix

    Const.use_CPP_engine = false; % 2018.07.31: Added now a flag to enable fast C++ Z matrix engine calculator
                            
    % Note:other parameters like Omega, K, Lambda, is dependent on the wavelenght and 
    % has to be calculated inside the frequency loops associated with the various routines, 
    % e.g. impedance matrix filling and is therefore not set here.

    % =================================
    % Some file-format version numbers that we support
    % =================================    
    Const.FEKO_efe_file_format = 4; % FEKO *.efe file format
    % TO-DO: Also check other files, e.g. *.str, *.hfe, *.ffe, *.snp, etc.    

    % Check now whether domain decomposition is active. This is signalled 
    % by the type of solver used, e.g. DGFM, CBFM, etc.
    if (Const.runCBFMsolver || Const.runJacobisolver || Const.runIFBMoMsolver ||Const.runDGFMsolver)
        message_fc(Const,sprintf('  Domain decomposition active'));
        Const.domain_decomposition = true;
    else
        Const.domain_decomposition = false;
    end

    % Finished
    message_fc(Const,sprintf('Done'));