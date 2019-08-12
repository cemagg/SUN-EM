function writeCMAsolToFile (Const, cma, xVectors)
    %writeCMAsolToFile v0.1
    %   Date: 03.04.2014
    %   Usage:
    %       writesolToFile(Const, CMA)
    %
    %   Input Arguments:
    %       Const
    %           A global struct containg the parameters, e.g., debug
    %           settings, *.str filename, etc.
    %       cma
    %           The CMA solution structure
    %
    %   Output Arguments:
    %           --
    %
    %   Description:
    %       Reads in a FEKO binary X sol. vector / *.str file which contains 
    %       either one or more set of MoM expansion coefficients. The Xsol 
    %       calculated with the CMA solution is then written to this file, so 
    %       that it can be viewed in POSTFEKO.
    %
    %   =======================
    %   Written by Danie Ludick on April 03, 2014
    %   Last updated on April 03, 2014.
    %   EM Systems & Software S.A. (Pty) Ltd.
    %   Email: dludick@emss.co.za
    
    %   References:
    %   [1] http://www.feko.info/support/helpcenter/how-to/how-to-read-the-.mat-.
    %       lud-.rhs-files-and-.str-files
    
    
    %   General hints for the reparseable format (option -r, i.e. FEKO can import again):
    %          - We write an ASCII *.str file in version 4.
    %          - The structure of this file is as follows:
    %                 4                                  <= version
    %                ################################    <= MD5
    %                n_MoM                               <= number of MoM currents
    %                n_FEM                               <= number of FEM currents
    %                --- end of header ---               <= header separation line
    %                (-###,-###)
    %                ...                                 <= n_MoM data lines
    %                (-###,-###)
    %                (-###,-###)
    %                ...                                 <= n_FEM data lines
    %                (-###,-###)
    %                --- separation of different blocks ---
    %
    %         - If there are several solutions, then the block with the data lines
    %           and the sparation line is repeated multiple times.
    %         - The separation line is only printed to separate the different solutions.
    %           (No separator between MoM and FEM currents!)
        
    error(nargchk(3,3,nargin));
        
    filename = Const.FEKOcmastrfilename;

    % Set the inpu filename
    fidIn  = fopen(sprintf('ascii_%s',filename),'r');
    
    % Set the output filename - depends on where this is called from (DGFM or standard CMA)
    if (Const.runCMAfromDGFMsolver || Const.useZactCMA)
        % We are running the CMA solver from the DGFM - append here also the array
        % element number
        % Overwrite the current *.str file
        fidOut = fopen(sprintf('cma_out_ascii_%s_array_element_%d',filename,Const.arrayID),'w');
    else
        % Normal CMA solver run - just use the filename as specified
        fidOut = fopen(sprintf('cma_out_ascii_%s',filename),'w');
    end%if    
    % Set the *.str filename

    if ((fidIn == -1) || (fidOut == -1))
        error(['Error writing the *.str file for FEKO']);
        mesage(Const,sprintf('Error writing the *.str file for FEKO'));
    end
    
    message_fc(Const,sprintf('Writing FEKCMA %s solution to *.str file: %s',upper(cma.name), filename));
    
    % Read and write the header information from the input file to the output file

    % Fileversion
    tline = fgets(fidIn);
    fprintf(fidOut,tline);
    
    % MD5 Checksum
    tline = fgets(fidIn);
    fprintf(fidOut,tline);
    
    % Number of MoM expansion coefficients (global)
    tline = fgets(fidIn);
    numMoMbasis=sscanf(tline,'%d');
    % Check that this is consistent with the number that was calculated by the CMA solver
    if (numMoMbasis ~= size(cma.Isol,1))
        error(['Error writing %s *.str file for FEKO - Isol length discrepency' upper(cma.name)]);
        mesage(Const,sprintf('Error writing %s *.str file for FEKO - Isol length discrepency', upper(cma.name)));
    end%if    
    fprintf(fidOut,tline);
   
    % Number of FEM expansion coefficients - not used
    tline = fgets(fidIn);
    fprintf(fidOut,tline);
    
    % The end-of-header delimiter
    tline = fgets(fidIn);
    fprintf(fidOut,tline);
    
    % Now write the eigencurrents to the *.str file
    %totalSols = cma.numSols*cma.numFreq;
    
    for freq=1:cma.numFreq
    	for ii=1:Const.numCMAmodes
            for jj=1:numMoMbasis
                fprintf ( fidOut, '(%20.10e,%20.10e)\n', real(cma.Isol(jj,ii,freq)), imag(cma.Isol(jj,ii,freq)) );
            end
            % write a delimiter line that seperates the individual modes    
            fprintf(fidOut,'--- separation of different blocks ---\n');
        end
    end
    
    fclose(fidIn);
    % delet the temporary input file
    %delete(sprintf('ascii_%s',filename));
    fclose(fidOut);    
    
    message_fc(Const,sprintf('  Wrote: %d solutions (number of Xsols.)',cma.numFreq*Const.numCMAmodes));
    %message_fc(Const,sprintf('  Wrote: %d solutions (number of Xsols.)',cma.numSols));
    message_fc(Const,sprintf('  Finished processing the *.str file: %s',filename));