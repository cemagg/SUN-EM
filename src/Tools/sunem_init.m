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
        if (Const.runCBFMsolver)
            % It is set in the driver
            set = true;
        end
    catch
        if (~set)
            % No CBFM solver present
           Const.runCBFMsolver = false;
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

    % =================================
    % Other initialisations:
    Const.useACA = false; % ACA interface available, but not yet activated. (Requires some refactoring).
    Const.ACAalg = 3; % When ACA is activated, we use the efficient algorithm (3).
    Const.store_to_fcd_file = false; % Store to external *.fcd file - useful when importing into POSTFEKO.
    Const.useCSCBFM = false;
    Const.isPhasedArray = false;

    % =================================
    % General constants:
    Const.RAD2DEG = 180/pi;  % For converting between radians and degrees
    Const.DEG2RAD = pi/180;  % For converting between degrees and radians

    Const.ETA = 120*pi;      % Wave-impedance of freespace
    Const.C0  = 3*10^8;      % Speed of light in free space
    Const.EPS = 10e-6;       % Set the machine precision to be used

    % =================================
    % Some file-format version numbers that we support
    % =================================    
    Const.FEKO_efe_file_format = 4; % FEKO *.efe file format
    % TO-DO: Also check other files, e.g. *.str, *.hfe, *.ffe, *.snp, etc.    

    % Finished
    message_fc(Const,sprintf('Done'));