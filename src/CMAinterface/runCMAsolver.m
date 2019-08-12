function [cma] = runCMAsolver(Const, zMatrices, yVectors)
    %runCMAsolver v0.1
    %   Date: 03.04.2014
    %   Usage:
    %       [cmaStruct] = runCMAsolver(Const, zMatrices)
    %
    %   Input Arguments:
    %       Const   
    %           A global struct containg the parameters, e.g., debug
    %           settings, *.str filename, etc.
    %       zMatrices
    %           The info related to the Z matrix as read from the binary
    %           FEKO file (e.g. M, N, precision, etc.)
    %       yVectors
    %           The Yrhs-vector data    
    %
    %   Output Arguments:
    %       cmaStruct
    %           Struct containing information such as the characteristic 
    %           currents and the eigenvalues. 
    %           Note: The Characteristic modes [J1 J2 .. Jn] for n = numModes 
    %           are sorted from smallest to largest absolute eigenvalue
    %
    %   Description:
    %       Performs a characteristic mode analysis on the complex
    %       impedance matrix, Z and returns the (sorted) eigenvalues and
    %       eigenvectors
    %
    %   =======================
    %   Written by Danie Ludick on March 26, 2013
    %   Last updated on April 04, 2014.
    %   EM Software & Systems - S.A. (Pty) Ltd.
    %   Email: dludick@emss.co.za
    
    %   References:
    %   ----------
    %   [1] Ludick, D. J., E. Lezar, and U. Jakobus. "Characteristic mode analysis of arbitrary 
    %       electromagnetic structures using FEKO." Electromagnetics in Advanced Applications 
    %       (ICEAA), 2012 International Conference on. IEEE, 2012.
    
    error(nargchk(3,3,nargin));    
    
    numSols = yVectors.numRhs;            % The number of right hand sides (excitations)
    % We need to run the DGFM (or i-DGFM) for each frequency point.
    numFreq = zMatrices.numFreq;           % The number of frequency points to 
    numRHSperFreq = numSols / zMatrices.numFreq;
                                           % The number of solutions per frequency point
    
    if (~Const.runCMAfromDGFMsolver)
        message_fc(Const,' ');
        message_fc(Const,'------------------------------------------------------------------------------------');
        if (Const.useZactCMA)
            message_fc(Const,sprintf('Running CMA-DGFM solver'));
            message_fc(Const,sprintf('Processing array element: %d',Const.arrayID));
            
            % Some info about the solution configurations
            message_fc(Const,sprintf('  numSols : %d', numSols));
            message_fc(Const,sprintf('  numFreq : %d', numFreq));
            message_fc(Const,sprintf('  numRHSperFreq : %d', numRHSperFreq));
            
            % Check now that the number of RHS per frequency point is 1. At
            % the moment, we only allow this for the CMA-DGFM
            
        else
            message_fc(Const,sprintf('Running CMA solver'));
        end
        % Write out some additional information:
        message_fc(Const,sprintf('Number of modes requested: %d',Const.numCMAmodes));
        message_fc(Const,sprintf('Number of frequencies    : %d',zMatrices.numFreq));
    else
        % Running the CMA from the DGFM solver
        message_fc(Const,' ');
        message_fc(Const,'    ------------------------------------------------------------------------------------');
        message_fc(Const,sprintf('  Running CMA solver from DGFM'));
        % Write out some additional information:
        message_fc(Const,sprintf('    Number of modes requested: %d',Const.numCMAmodes));
        message_fc(Const,sprintf('    Number of frequencies    : %d',zMatrices.numFreq));
        message_fc(Const,sprintf('    Array element ID         : %d',Const.arrayID));
    end%if

    % Start the CMA timing
    % TO-DO: What to do when this is called from within the DGFM?
    tic        
    
    % Extract info related to the Z matrix    
    M = zMatrices.mBasis;
    N = zMatrices.nBasis;
    numFreq = zMatrices.numFreq;
    %numModesRequested = Const.numCMAmodes;   % Not required here anymore (see below)
    if (Const.checkRXorth) 
        cma.numOrthModes = zeros(numFreq,1);
    end%if

% The following is now done inside the driver routine, before runCMA is called
%     % Calculate the number of internal modes to calculate by adding a
%     % buffer percentage - similar to KERNEL (only if tracking is used)
%     if (Const.trackModes)
% %     -- Tracking used, set the percentage increase in the number of modes that will be calculated
% %        Note, if the used requested less than 10 modes, then double the number of modes requested
%         if (numModesRequested <= 5)
%         % -- less than 5 modes - use 100% increase
%           TRACKING_BUFF_PERCENTAGE = 1.0D0;
%         elseif ( (numModesRequested >= 5) && (numModesRequested <= 10) )
%         % -- between 5 and 10 modes - use 50% increase
%           TRACKING_BUFF_PERCENTAGE = 0.5D0;
%         else
%         % -- 10 modes and more - use 25% increase
%           TRACKING_BUFF_PERCENTAGE = 0.25D0;
%         end%if
%     else
%         TRACKING_BUFF_PERCENTAGE = 0.0D0;
%     end%if
    
    % Increase now the number of modes with the tracking % buffer
    %numModes = numModesRequested + ceil(numModesRequested*TRACKING_BUFF_PERCENTAGE);
    numModes = Const.numModesCalc;
    message_fc(Const,sprintf('Number of modes calculated: %d',numModes));

    % Loop over each of the frequencies
    for freq = 1:numFreq

        message_fc(Const,sprintf('  Processing frequency index: %d of %d',freq, numFreq));
        
        % Extract Z matrix at current frequency
        Z = zMatrices.values(:,:,freq);
	        
        % Symmetrise the Z-matrix if requested. Note: In FEKO, this is not done
        % explicitely, as numerous tests showed that it is not that critical.
        if (Const.symZmat)
            message_fc(Const,sprintf('    Calculating the symmetric matrix'));
            % Loop over the upper part of the matrix Z
            limit = 0.01D0;
            for jj = 1:M
                for ii = jj+1:N
                % calcalate the average of the entries Zij and Zji
                    diff = (Z(jj,ii) - Z(ii,jj));
                    ctemp = 0.5D0*( Z(jj,ii) + Z(ii,jj) );
                % calculate the relative difference in the matrix elements Zij and Zji    
                    temp  = abs(  diff / ctemp );
                % check whether the relative difference is less that our predefined limit
                    if ((temp > limit) && (Const.debug))
                % this is an error    
                        message_fc(Const,sprintf ('Warning: Z(%d,%d) and Z(%d,%d), Rel. diff of %6.2f',ii,jj,jj,ii,temp));
                    end
                % set both elements now equal to the average
                    Z(jj,ii) = ctemp;
                    Z(ii,jj) = ctemp;
                end
            end
        end
	
        % Test here whether the matrix Z is symmetrical
        if ( Const.symZmat && (~isequal(Z,Z.')) )
            mesage(Const,sprintf('Z is not a symmetrical matrix'));
            error(['Z is not a symmetrical matrix']);
        end
	
        % Calculate the real and imaginary parts of the impedance matrix
        R = -1.*real(Z);
        X = -1.*imag(Z);
	
        % Test here whether the matrix Z is symmetrical
        if ( Const.symZmat && ((~isequal(R,R.')) || (~isequal(X,X.'))) )
            mesage(Const,sprintf('R or X is not a symmetrical matrix'));
            error(['R or X is not a symmetrical matrix']);
        end 
	
        % ====================================================
        % ===   Solve the generalised eigenvalue equation: ===
        % ===           [X][I]n  = \lambda[R][I]n          ===
        % ====================================================
	
        % Using the eig-function. Use 'chol' factorisation, as we are working with 
        % symmetrical matrices X and R    
        
        if (Const.eigSolver == 1)
            % LAPACK routines
            message_fc(Const,sprintf('    Solving the eigenvalue problem using LAPACK'));
            [EigenVec, EigenValMat] = eig(X,R,'chol');
            
        elseif (Const.eigSolver == 2)
            message_fc(Const,sprintf('    Solving the eigenvalue problem using ARPACK (SM)'));        
            % ARPACK (SEP)
            % First convert the generalised eigenvalue problem to a standard
            % eigenvalue problem and then solve for the smallest eigenvalues
            D = R\X;
            [EigenVec, EigenValMat] = eigs(D, numModes, 'SM');
            
        elseif (Const.eigSolver == 3)
            message_fc(Const,sprintf('    Solving the eigenvalue problem using ARPACK (LM)'));                
            % ARPACK (SEP)
            % First convert the generalised eigenvalue problem to a standard
            % eigenvalue problem (other way around) and then solve for the 
            % largest eigenvalues (other way around, which is more stable
            % apparantly) - See [1].
            D = X\R;
            % For debug information we activate some additional ARPACK
            % output:
            if (~Const.debug)
                OPTS.disp = 0;
            else
                OPTS.disp = 2;
            end%if
            [EigenVec, EigenValMat] = eigs(D, numModes, 'LM', OPTS);
	
        elseif (Const.eigSolver == 4)
            % LAPACK (SEP)
            % Solve the GEP and find eigenvalues closest to sigma (in our case 0)
            sigma = 0; % Find eigenvalues close to sigma
            % In order for the GEP solver to work, [R] has to be
            % positive semi-definite. In theory this is the case, but due to
            % numerical tolerances this is seldomly the case. In [2],
            % Harrington discusses a method to force [R] to be positive
            % definite by doing an SVD and then omitting some of the smaller
            % and negative eigenvalues. This type of preconditioning is
            % investigated here.
            [U,EigenValMat,V] = svd(R);
            EigenValTmp = diag(EigenValMat,0);
            
            % Just perform an intermediary check here to make sure that the
            % u and v vectors are correct (i.e. can be used as in [2])
            EigenValCheck = (U.')*R*U;        
            diagCheck = diag(EigenValCheck,0);
            err = calculateErrorNormPercentage(EigenValTmp, diagCheck);
            if (err >= (10^(-3)) )
                message_fc(Const,sprintf('Error during R-preconditioning'));
                error(['Error during R-preconditioning']);
            end%if
            
            % Plot the eigenvalues on a log scale
            figure;
            hold on;
            grid on;
            box on;
            plot([1:length(EigenValTmp)], log10(EigenValTmp), 'r-*');
            threshPlot = [1:length(EigenValTmp)];
            threshPlot(:) = Const.RorthThresh;
            plot([1:length(EigenValTmp)],log10(threshPlot), 'k--');
            
            % See [2] Eq. (18): Partition diagonal matrix by removing small
            % eigenvalues:
            count = 0;
            A22_offset_set = false;
            for ii=1:length(EigenValTmp)
                if (EigenValTmp(ii) < Const.RorthThresh)
                    if (~A22_offset_set)
                        A22_offset = ii;
                        A22_offset_set = true;
                    end%if
                    count = count+1;
                    % Calculate [mu] - see Eq. (18) in [2]
                    EigenValMat(ii,ii) = 0.0d0;
                end%if
            end%for
            message_fc(Const,sprintf ('    R-orthogonality: Zerod %d out of %d entries',count,length(EigenValTmp)));
            
            % See Eq. (20) in [2]
            A = U.'*X*U;
            %   From the A, extract the submatrices A11, A12, A12t and A22
            A11 = A(1:A22_offset-1,1:A22_offset-1);
            A12 = A(1:A22_offset-1,A22_offset:N);
            A12t = A(A22_offset:M,1:A22_offset-1);
            A22 = A(A22_offset:M,A22_offset:N);
            
            % TO-DO: Check whether A12 is in-fact A12t
            
            % Calculate now the real symmetric unweighted eigenvalue equation
            % (see Eqn. (26) in [2]), i.e., [B][y] = lambda[y] with B as follows:
            mu_11 = sqrt(EigenValMat(1:A22_offset-1,1:A22_offset-1));
            % Loop over the diagonal entries and calculate the inverse
            for ii=1:A22_offset-1
                mu_11_inv(ii,ii) = 1./mu_11(ii,ii);
            end%for
            B = mu_11_inv*(A11 - A12*inv(A22)*A12t)*(mu_11_inv);
            
            % Check whether we need to reduce now the number of modes that were requested:
            if (numModes >= (A22_offset-1))
                numModes = A22_offset-1;
                message_fc(Const,sprintf ('    R-orthogonality: System dimensions reduced to %d entries',numModes));
                
                % In this case we use LAPACK
                message_fc(Const,sprintf('    Solving the eigenvalue problem using LAPACK'));
                [EigenVec, EigenValMat] = eig(B);
                
            else
                % Calculate the number of modes that have been requested using ARPACK (SEP)
                message_fc(Const,sprintf('    Solving the eigenvalue problem using ARPACK (SM)'));
                [EigenVec, EigenValMat] = eigs(B,numModes,'SM');
	
            end%if
            
            % Calculate now the eigenvectors according to Eq. (27) in [2]
            tmp = -inv(A22)*A12t;
            delta = eye(A22_offset-1,size(tmp,2));
            EigenVec = U*[delta ; tmp ]*mu_11_inv*EigenVec;
            
        else
            message_fc(Const,sprintf('Unsupported CMA algorithm'));
            error(['Unsupported CMA algorithm']);
        end    
	
        % Extract now only the diagonals of the matrix EigenValMat
        EigenVal = diag(EigenValMat,0);
        if (Const.eigSolver == 3)
            % Invert the eigenvalues
            EigenVal(:) = 1./EigenVal; % Invert the eigenvalues
        end
	
        % We are interested in the eigenvectors (modes) with the lowest eigen-value
        % magnitude. Depending on the eigensolver used, we need to perform a
        % sorting here.
	
        % Print some matrix entries for debugging purposes (ony when debugging is active)
        if (Const.debug)    
            num_entries = 5;
            message_fc(Const,sprintf('    The first %d matrix entries:',num_entries*num_entries));
            % Print out [X]
            for ii=1:num_entries
                for jj =1:num_entries
                    message_fc(Const,sprintf('      AEigMat(%d,%d) = %.7E',ii,jj,X(ii,jj)));
                end
            end
            fprintf('\n');
            % Print out [R]
            for ii=1:num_entries
                for jj =1:num_entries
                    message_fc(Const,sprintf('      BEigMat(%d,%d) = %.7E',ii,jj,R(ii,jj)));
                end
            end
        end
        
        % Sort the eigenvalues if they were calculated with LAPACK
        % Always sort the eienvalues and eigenvectors - seeing the same
        % behaviour in MATLAB 6.5.0.180913a (R13) with ARPACK not sorting
        % the values correctly
        %if ( Const.eigSolver == 1 || ((Const.eigSolver == 4) && (numModes < (A22_offset-1))) )
            % Sort the Eigenvalues accordingly and also extract a permutation vector 
            [EigenValS, Perm] = sort(abs(EigenVal));
            EigenValS = EigenVal(Perm(:));
            EigenVal = EigenValS;
            % Use now the permutation vector to sort also the eigen-currents that is stored
            % in the matrix EigenVec
            EigenVecS = EigenVec(:,Perm(:));
            EigenVec  = EigenVecS;
        %end
            
        % Print the number of requested eigen-values to the screen (only in debug mode)
        if (true) 
            message_fc(Const,sprintf('      The first %d eigen-modes:',numModes));
            for ii=1:numModes
                message_fc(Const,sprintf('        mode-indes: %d mode-val: %.5E',ii,EigenVal(ii)));
            end
        end
        
        % Normalise modes according to the maximum value
        for mm=1:numModes        
            MaxValue = max(EigenVec(:,mm));
            EigenVec(:,mm) = EigenVec(:,mm)./MaxValue;
        end
        
        % Normalise the EivenVectors such that <Jn*,RJn> = 1
        for mm=1:numModes
            norm = (EigenVec(:,mm).')*(R*EigenVec(:,mm));
            sign = +1;
            if (norm >= 0) 
                sign = +1;
            else
                sign = -1;
            end
            % take the ABS of norm, otherwise sqrt(norm) will be imaginary
            norm = abs(norm);
            EigenVec(:,mm) = sign*EigenVec(:,mm)./sqrt(norm);
        end
        
        % Perform some checks for R and X-orthogonality of the modes
        if (Const.checkRXorth) 
            cma.numOrthModes(freq,1) = checkModeOrthogonality(Const, EigenVal, EigenVec, R, X);
        end%if
        
        cma.solTime = toc;
        
        % Set the return values
        cma.name    = 'CMA';
        cma.Isol(:,:,freq) = EigenVec;
        cma.EigVal(:,freq)  = EigenVal;
        %cma.numSols = Const.numCMAmodes;
        cma.numSols = numModes;
        cma.numFreq = zMatrices.numFreq;
%       cma.Rorth = Rorth;
%       cma.Xorth = Xorth;
    
%       See issue FEKDEV-25040: Perform tracking (calculate a mod_rank
%       to map the current eigenvalues to those at the previous frequency
%       increment) - if we follow the KERNEL tracking approach. Otherwise
%       for adaptive tracking we can easily do this after the frequency
%       increment.
%       TO-DO: I think adaptive tracking enabled us to track the modes AFTER 
%       the frequency loop. Therefore comment out that option below.
        if (Const.trackModes == 1)
            [cma] = cmaTrackModes (Const, cma, EigenVal, EigenVec, R, X, freq);
        %elseif (Const.trackModes == 2)
        %    [cma] = cmaTrackModesAdaptive(Const, cma);
        end%if
    
    end%for freq = 1:numFreq
    
    % Also calculate and set the characteristic angle:
    cma.CharAngle = zeros(size(cma.EigVal));
    for mm = 1:numModes
        cma.CharAngle(mm,:) = 180 - atan(cma.EigVal(mm,:))*180/pi;
    end%for
    
    message_fc(Const,sprintf('Finished the CMA solver in %f sec.',cma.solTime));
