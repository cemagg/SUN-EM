function writeSolToFile(Const, solStruct)
    %writeCBFMsolToFile v0.1
    %   Date: 25.06.2013
    %   Usage:
    %       writesolToFile(Const, cbfm)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containging:
    %       filename
    %           FEKO *.str filename (e.g. 'cbfm.str') (as output by FEKO)
    %           The expansion coefficients in this file is replaced with
    %           those calculated by FEKDDM's DD solution (e.g. the CBFM,
    %           DGFM, Jacobi Iterative solver, Iterative DGFM, etc.).
    %       solStruct
    %           Structs containing details about the solver - and which
    %           *.str filename to use.
    %
    %   Output Arguments:
    %           --
    %
    %   Description:
    %       Reads in a FEKO binary X sol. vector / *.str file which contains
    %       either one or more set of MoM expansion coefficients, each at a
    %       different frequency / solution configuration. The Xsol calculated
    %       with FEKDDMs DD solution (e.g. CBFM) is then written to this file, so
    %       that it can be viewed in POSTFEKO.
    %
    %   =======================
    %   Written by Danie Ludick on June 25, 2013
    %   Last updated on June 19, 2011.
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

    error(nargchk(2,2,nargin));

    % Make sure to use the correct filename
    no_output = false;
    switch (upper(solStruct.name))
        case 'CBFM'
            filename = Const.SUNEMcbfmstrfilename;
        case 'DGFM'
            filename = Const.SUNEMdgfmstrfilename;
        case 'JACK'
            filename = Const.SUNEMjackstrfilename;
        case 'NGFDGFM'
            filename = Const.SUNEMngfdgfmstrfilename;
        case 'IFBMOM'
            filename = Const.SUNEMifbmomstrfilename;
        case 'ITERDGFM'
            filename = Const.SUNEMiterdgfmstrfilename;
        case 'HARP-MBF-PATTERN'
            filename = Const.harp_single_element_MBF_patterns_FEKOstrfile;
            no_output = true;
        otherwise
            message(Const,sprintf('Unsupported filename: ',solStruct.name));
            error(['Unsupported filename']);
    end

    ascii_strfilename = sprintf('ascii_%s',filename);
    if (exist(ascii_strfilename,'file'))
        % Do nothing the ASCII file will be read later
    else
    % First convert the *.str file to a ASCII file (only on Windows, see FEKDDM-08)
        if (ispc)
            [stat, dosout] = dos(sprintf('\"%s//str2ascii.exe\" %s -r > %s 2>&1', ...
                Const.ExecPath, filename, ascii_strfilename));
            if (stat ~= 0)                
                message(Const,sprintf('Error converting FEK0 *.str file to ASCII: %s',filename));
                message(Const, num2str(stat));
                error(['Error converting FEKO *.str file to ASCII: %s' filename]);
            end%if

        elseif (isunix)
            [stat, sysout] = system(sprintf('\"%s//str2ascii\" %s -r > %s 2>&1', ...
                Const.ExecPath, filename, ascii_strfilename));
            %sysout
            if (stat ~= 0)
                sysout % Write the FEKO output to the screen
                error(['Error running str2ascii for: ' filename]);
                message_fc(Const,sprintf('Error running str2ascii for: %s', filename));
                message_fc(Const, num2str(stat));
            end%if

        else
		% See FEKDDM-08, the str2ascii utility is not working on MAC iOS. We need
	    % to work there with an already converted *.str file			
			message_fc(Const,sprintf('Error: Converting FEKO *.str file to ASCII only possible on Windows'));
            error(['Error: Converting FEKO *.str file to ASCII only possible on Windows']);
	    end%if
    end%if (exist(ascii_strfilename,'file'))

    fidIn  = fopen(sprintf('ascii_%s',filename),'r');
    % Overwrite the current *.str file
    fidOut = fopen(sprintf('out_ascii_%s',filename),'w');
    if ((fidIn == -1) || (fidOut == -1))        
        message_fc(Const,sprintf('Error writing the *.str file for FEKO'));
        error(['Error writing the *.str file for FEKO']);
    end

    if (~no_output)
        message_fc(Const,sprintf('  Writing %s solution to *.str file: %s',upper(solStruct.name), filename));
    end%if

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
    % Check that this is consistent with the number that was calculated by the FEKODDM DD solver (e.g. the CBFM or the DGFM)
    if (numMoMbasis ~= length(solStruct.Isol))
        error(['Error writing %s *.str file for FEKO - Isol length discrepency' upper(solStruct.name)]);
        message(Const,sprintf('Error writing %s *.str file for FEKO - Isol length discrepency', upper(solStruct.name)));
    end%if
    fprintf(fidOut,tline);

    % Number of FEM expansion coefficients - not used
    tline = fgets(fidIn);
    fprintf(fidOut,tline);

    % The end-of-header delimiter
    tline = fgets(fidIn);
    fprintf(fidOut,tline);

    % Now write the DD expansion functions to file (for each solution configuration)
    numSols = solStruct.numSols;
	for solNum=1:numSols
        for jj=1:numMoMbasis
            fprintf ( fidOut, '(%20.10e,%20.10e)\n', real(solStruct.Isol(jj,solNum)), imag(solStruct.Isol(jj,solNum)) );
        end
	    % write a delimiter line that seperates the individual modes
        fprintf(fidOut,'--- separation of different blocks ---\n');
    end%for

    fclose(fidIn);
    % Do not delete the temporary input file, as the str2ascii utility
    % cannot be used on certain platforms, such as Mac OS X.
    %delete(sprintf('ascii_%s',filename));
    fclose(fidOut);

    if (~no_output)
        message_fc(Const,sprintf('  Wrote: %d solutions (number of Xsols.)',numSols));
        message_fc(Const,sprintf('  Finished processing the *.str file: %s',filename));
    end%if