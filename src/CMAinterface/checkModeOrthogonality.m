function [numOrthModes] = checkModeOrthogonality(Const, EigenVal, EigenVec, R, X)
    %runCMAsolver v0.1
    %   Date: 26.05.2014
    %   Usage:
    %       [cma] = checkModeOrthogonality(Const, EigenVal, EigenVec, R, X)
    %
    %   Input Arguments:
    %       Const   
    %           A global struct containg the parameters, e.g., debug
    %           settings, *.str filename, etc.
    %       EigenVal, EigenVec
    %           The eigen values and eigenvectors as calculated in
    %           runCMAsolver.m
    %       R, X
    %           The R and X matrices (real and imaginary components of the
    %           MoM impendance matrix such that Z = R + j*X)
    %
    %   Output Arguments:
    %       numOrthModes
    %           The number of R and X orthogonal modes
    %
    %   Description:
    %       Calculates R and X orthogonality and prints out additional
    %       information if requested (see Const.checkRXorth):
    %       Const.checkRXorth  = 0: Nothing
    %                            1: Calculate the number of orthogonal
    %                               modes (same as in KERNEL)
    %                            2: Do (1) and write out the orthogonality 
    %                               of each mode
    %                            3: Do (1) and (2) and plot the R and X
    %                               orthogonality                    
    %
    %   =======================
    %   Written by Danie Ludick on May 26, 2014
    %   Last updated on May 26, 2014.
    %   EM Software & Systems - S.A. (Pty) Ltd.
    %   Email: dludick@emss.co.za
    
    %   References:
    %   ----------
    
    error(nargchk(5,5,nargin));

    message(Const,sprintf('  Perform R and X orthogonality check'));

    % Initialise variables
    CMA_MODE_ORTH_THRESH = 0.01D0;
    numModes = Const.numCMAmodes;
    numOrthModes = 0;
    
    % Perform same orthogonality check as in FEKO KERNEL
    if (Const.checkRXorth >= 1)
        
%     Initialise here the number of modes that can be used to 
%     the maximum amount
        numcmamodes_limit = numModes;
%     Loop over all the modes and calculate <Jm, [R]Jn> and <Jm, [X]Jn>
        for nn=1:numModes
%           Get the maximum number of orthogonal modes
            if (nn > numcmamodes_limit)
%             Distiguish between two cases here to obtain the maximum number of modes
              if ((nn-1)==numcmamodes_limit)
%               Special case here where we should limit the number of modes to JN_IND-1
                numcmamodes_limit = nn - 1;
              else
%               Limit the modes to the maximum number according to the N-index. Note, we subtract 2
%               here and not 1 (as done for the JM_IND case below), as JN_IND is incremented above in the
%               CMA_JN_LOOP before we get here.
                numcmamodes_limit = nn - 2;
              end%if
%             Exit from checks - return this value as the number of orthogonal modes
              %EXIT CMA_JN_LOOP
              break;
            end%if
            
%           Extract the eigencurrent, Jn            
            Jn = EigenVec(:,nn);
            RJn = (R*Jn);
            XJn = (X*Jn);
            
%           Compare now the product [R]*Jn with each of the modes, Jm            
            for mm=1:numModes
              if (mm > numcmamodes_limit )
%               Exit from loop if index exceeds our limit
                break;
              end%if
%             First extract Jm              
              Jm = EigenVec(:,mm);
              
%             Calculate the mode orthogonalities
              norm_r = (Jm.')*RJn;
              norm_x = (Jm.')*XJn;
              
              if (mm == nn)              
%             -- <Jm,[R]Jn> = 1 if m = n and <Jm,[X]Jn> = EigVal if m = n
%                otherwise the modes are not orthogonal
                if ( (abs(1.0d0 - norm_r) >= CMA_MODE_ORTH_THRESH) || ...
                     (abs(real(EigenVal(mm)) - norm_x) >= CMA_MODE_ORTH_THRESH) )
%                   Set the limit according to the M index, as we cal still check up and to 
%                   <Jm,(R/X)Jm> i.e. N=M.                 
                    numcmamodes_limit = mm - 1;
                    message(Const, sprintf('--- Limiting mode Jm to: %d',numcmamodes_limit));
                end%if

              else
%             -- <Jm,[R]Jn> = 0 if m <> n and <Jm,[X]Jn> = 0 if m <> n

                if ( (abs(norm_r) >= CMA_MODE_ORTH_THRESH) || ...
                     (abs(norm_x) >= CMA_MODE_ORTH_THRESH) )
%                   Set the limit according to the M index (see above comment why we do this)                
                    numcmamodes_limit = mm - 1;
                    message(Const, sprintf('--- Limiting mode Jm to: %d',numcmamodes_limit));
                end%if
              end%if
              
            end%for mm=1:numModes
        end%for nn=1:numModes

% If the limit here is zero, then we only consider the first mode for the remainder of the 
% calculations
        if (numcmamodes_limit < 0) 
            numcmamodes_limit = 1;
        end%if

% Set the return value
        numOrthModes = numcmamodes_limit;        
        message(Const, sprintf('Number of orthogonal eigenmodes: %d',numOrthModes));
        
    end%if (Const.checkRXorth >= 1)
    
    % Write out R and X-orthogonality of the modes
    if (Const.checkRXorth >= 2)
%        Rorth = zeros(numModes,numModes);
%        Xorth = zeros(numModes,numModes);
        for nn=1:numModes
            Jn = EigenVec(:,nn);
            RJn = (R*Jn);
            XJn = (X*Jn);
            for mm=1:numModes
                Jm = EigenVec(:,mm);            
                message(Const,sprintf('    <J%d,RJ%d> = %.5E',mm,nn,(Jm.')*RJn));
                message(Const,sprintf('    <J%d,XJ%d> = %.5E',mm,nn,(Jm.')*XJn));
%                Rorth(mm,nn) = (Jm.')*RJn;
%                Xorth(mm,nn) = (Jm.')*XJn;
            end
        end        
    end%if
        
    if (Const.checkRXorth >= 3)
        % Plot the R-orthogonality of modes (set here which mode index nn)
        nn = 1;
        figure;
        semilogy([1:numModes],abs(Rorth(1:numModes,nn)));
        title(sprintf('R-orthogonality of <Jm,RJ%d> with m=%d to %d',nn,1,numModes));
        xlabel('Mode index (m)');
        ylabel('10log(<Jm,RJn>)');
        grid on;
        box on;
        
        % Plot the X-orthogonality of modes
        nn = 1;
        figure;
        semilogy([1:numModes],abs(Xorth(1:numModes,nn)));
        title(sprintf('X-orthogonality of <Jm,XJ%d> with m=%d to %d',nn,1,numModes));
        xlabel('Mode index (m)');
        ylabel('<Jm,XJn>');
        grid on;
        box on;
    end%if
