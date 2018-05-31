function [xVectors] = readFEKOXvectorFromFile(Const, strfilename)
    %readFEKOZMatrixFromFile
    %   Usage:
    %       [xVectors] = readFEKOXvectorFromFile(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct containing program flow parameters and constants
    %       strfilename
    %           FEKO str filename (e.g. 'yagi.str') (as output by FEKO)
    %
    %   Output Arguments:
    %       xVectors, a struct containing:
    %       values
    %           Set of P FEKO X vectors (MxP), representing an M-dim. vector
    %           at P frequencies / solution configurations.
    %       fileVersion
    %           The file version of the *.str file
    %       md5Check
    %           The MD5 checksum of the *.str file
    %       numMoMbasis
    %           The number of MoM basis-functions
    %       numFEMbasis
    %           The number of FEM basis-functions
    %       numSols
    %           The number of solutions vectors
    %   Description:
    %       Reads in a FEKO binary X sol. vector / *.str file which contains
    %       either one or more set of MoM expansion coefficients, each at a
    %       different frequency / solution configuration and returns a
    %       struct with the relevant data.
    %
    %   References:
    %   [1] http://www.feko.info/support/helpcenter/how-to/how-to-read-the-.mat-.
    %       lud-.rhs-files-and-.str-files
    %
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

    % *.str file that is converted to ASCII using the str2ascii utility (see below).
    % If we already have such a file in the directory, then use it.
    ascii_strfilename = sprintf('ascii_%s',strfilename);

    if (exist(ascii_strfilename,'file'))
        % If we already have a *.str file that is converted to ascii - then use that.
    else
        % First convert the *.str file to a ASCII file (only on Windows and Linux)
        if (ispc)
            [stat, dosout] = dos(sprintf('\"%s//str2ascii.exe\" %s -r > %s 2>&1', ...
                Const.ExecPath, strfilename, ascii_strfilename));
            if (stat ~= 0)                
                message_fc(Const,sprintf('Error converting FEK0 *.str file to ASCII: %s', ...
                    strfilename));
                message_fc(Const, num2str(stat));
                error(['Error converting FEKO *.str file to ASCII: %s' strfilename]);
            end%if
        elseif (isunix)
            [stat, sysout] = system(sprintf('\"%s//str2ascii\" %s -r > %s 2>&1', ...
                Const.ExecPath, strfilename, ascii_strfilename));
            if (stat ~= 0)                
                message_fc(Const,sprintf('Error converting FEK0 *.str file to ASCII: %s', ...
                    Const.FEKOstrfilename));
                message_fc(Const, num2str(stat));
                error(['Error converting FEKO *.str file to ASCII: %s' strfilename]);
            end%if
        else
    		% See FEKDDM-08, the str2ascii utility is not working on MAC iOS. We need
	        % to work there with an already converted *.str file			
			message_fc(Const, ...
                sprintf('Error: Converting FEKO *.str file to ASCII only possible on Windows or Linux'));
            error(['Error: Converting FEKO *.str file to ASCII only possible on Windows or Linux']);
	    end%if
    end%if (exist(ascii_strfilename,'file'))


    % Open now this file
    fid = fopen(ascii_strfilename,'r');
    if fid == -1
        error(['Error reading FEKO *.str file: ' ascii_strfilename]);
        message_fc(Const,sprintf('Error reading: ascii_%s',ascii_strfilename));
    end

    message_fc(Const,' ');
    message_fc(Const, ...
        '------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Reading *.str file from: %s',ascii_strfilename));

    % Read the header information from the input file and write this to the
    tline = fgets(fid);
    xVectors.fileVersion = sscanf(tline,'%f');

    xVectors.md5Check = fgets(fid);

    tline = fgets(fid);
    xVectors.numMoMbasis=sscanf(tline,'%d');

    tline = fgets(fid);
    xVectors.numFEMbasis=sscanf(tline,'%d');

    % The end-of-header delimiter
    tline = fgets(fid);

    % Now follow the numMoMbasis basis-functions
    xVectors.numSols = 0;
    xVectors.values = [];
    while ~feof(fid)
        xVectors.numSols = xVectors.numSols + 1;
        for jj=1:xVectors.numMoMbasis
            tline = fgets(fid);
            tmp = sscanf(tline,'(%f,%f)');
            xVectors.values(jj,xVectors.numSols) = tmp(1) + 1i*tmp(2);
        end%for
        % read the delimiter line that seperates the individual modes
        tline = fgets(fid);
    end%for

    % Close the file again
    fclose(fid);

    message_fc(Const,sprintf('Read Xsol file version %d',xVectors.fileVersion));
    message_fc(Const,sprintf('No. of MoM basis = %d and no. of FEM basis = %d', ...
        xVectors.numMoMbasis,xVectors.numFEMbasis));
    message_fc(Const,sprintf('Read: %d solutions (number of Xsols.)',xVectors.numSols));
    message_fc(Const,sprintf('Finished processing the *.str file: %s',ascii_strfilename));