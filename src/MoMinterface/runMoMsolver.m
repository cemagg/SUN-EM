function [mom] = runMoMsolver(Const, Solver_setup, zMatrices, yVectors, refIsol)
    %runMoMsolver
    %   Usage:
    %       [mom] = runMoMsolver(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       Solver_setup
    %           Solver specific struct, e.g. frequency range, basis function details, geometry details    
    %       zMatrices
    %           The Z-matrices data. This can be from FEKO (extracted from the *.mat file, or internally
    %           calculated).
    %       yVectors
    %           The Yrhs-vector data
    %       refIsol
    %           The reference solution-vector data (e.g. MoM solution of FEKO or SUN-EM)
    %
    %   Output Arguments:
    %       mom
    %           Structs containing MoM solution and timing data
    %
    %   Description:
    %       Runs the MoM solution based on the Z and Y data that was read / parsed.
    %       from the FEKO *.out, *.mat, *.str and *.rhs files
    %
    %   =======================
    %   Written by Danie Ludick on June 21, 2013.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    narginchk(5,5);

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Running MoM solver'));

    % Initialise the return values
    mom  = [];
    mom.name = 'mom';
    Nmom = Solver_setup.num_mom_basis_functions;   % Total number of basis functions for whole problem
    numSols = refIsol.numSols;                     % The number of solutions configurations
    mom.numSols = numSols; %numSols;               % For now, set to 1. (TO-DO: Update)
    numFreq = Solver_setup.frequencies.freq_num;   % The number of frequency points to process
    numRHSperFreq = 1; % mom.numSols / numFreq;    % The number of solutions per frequency point.
                                                   % For now, should be 1 (TO-DO: Update)
                                                   
    % Some info about the solution configurations
    message_fc(Const,sprintf('  numSols : %d', mom.numSols));
    message_fc(Const,sprintf('  numFreq : %d', numFreq));
    message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));

    % Calculate the solution vector (all frequency points, all RHSes)
    mom.Isol = complex(zeros(Nmom,mom.numSols));

    % The timing calculations also need to take into account that there is a
    % frequency loop
    mom.setupTime = zeros(1,numFreq);
    mom.factorisationTime = zeros(1,numFreq);
    % Zero also the total times (for all frequency iterations)
    mom.totsetupTime = 0.0;
    mom.totfactorisationTime = 0.0;
    mom.totsolTime = 0.0;

    % Start the frequency loop now
    for freq=1:numFreq

        % Start timing (MoM setup)
        tic

        % Reset each solution per frequency point here
        solStart = 1;
        solEnd   = numRHSperFreq;

        % Allocate here space for the MoM matrix so that it can be filled
        ObservRWGs = [1:Nmom];
        SourceRWGs = [1:Nmom];
        % Note: Since 2017-06-25, we are also passing a freq. variable here to
        % indicate for which frequency point we are extracting the matrix
        Zmom = (calcZmn(Const, zMatrices, freq, 1, 1, ObservRWGs, SourceRWGs));

        % End timing (calculating the impedance matrix)
        mom.setupTime(freq) = toc;

        % Start timing (MoM factorisation)
        tic
        % LU-decomposition of the Z-matrix
        [L,U] = lu(Zmom);

        % Loop over the RHs vectors (numSols) and calculate each of the currents.
        for solNum=solStart:solEnd

            % Take care where we are extracting the values from (out of Yrhs)
            % and also where we will be storing these values again (in Xsol)
            index = solNum + (freq-1)*numRHSperFreq;
            % DJdbg --> remove
            %message_fc(Const,sprintf('  index: %d', index));

            % Back-wards substitution
            b = L\yVectors.values(:,index);
            mom.Isol(:,index) = U\b;
        end%for

        % End timing (MoM factorisation)
        mom.factorisationTime(freq) = toc;

        % Total time (MoM matrix setup + factorisation)
        mom.solTime(freq) = mom.setupTime(freq) + mom.factorisationTime(freq);

        % Memory usage remains constant between frequency iterations
        mom.memUsage = byteSize(Zmom);

        % Calculate the total MoM times
        mom.totsetupTime = mom.totsetupTime + mom.setupTime(freq);
        mom.totfactorisationTime = mom.totfactorisationTime + mom.factorisationTime(freq);
        mom.totsolTime = mom.totsolTime + mom.solTime(freq);

    end%for freq=1:numFreq

    message_fc(Const,sprintf('Finished MoM solver in %f sec.',mom.totsolTime));
    message_fc(Const,sprintf('(Times for Z-setup : %f sec. and LU-fact. : %f)',mom.totsetupTime, ...
        mom.totfactorisationTime));
    message_fc(Const,sprintf('Memory usage of MoM %s',mom.memUsage));

    % Compare the MoM solution obtained with MATLAB, with that obtained by FEKO
    % that was stored in xVectors.values (for each frequency iteration (and each solution within the frequency iteration)
    % Calculate also space for the relative error here
    mom.relError = zeros(1,mom.numSols);
    for freq=1:numFreq
        for solNum=1:numRHSperFreq
            index = solNum + (freq-1)*numRHSperFreq;
            mom.relError(index) = calculateErrorNormPercentage(refIsol.Isol(:,index), mom.Isol(:,index));
            message_fc(Const,sprintf('Rel. error norm. for Sol. %d of %d of freq. %d of %d compared to reference sol. %f percent',solNum, ...
                numRHSperFreq, freq, numFreq, mom.relError(index)));
        end
    end
    
    % Write the MoM solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested (i.e. if the filename is defined)
    if (~isempty(Const.SUNEMmomstrfilename))
        writeSolToFile(Const, mom);
    end%if

    