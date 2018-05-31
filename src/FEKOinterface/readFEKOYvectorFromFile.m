function [yVectors] = readFEKOYvectorFromFile(Const, rhsfilename)
    %readFEKOZMatrixFromFile
    %   Usage:
    %       [yVectors] = readFEKOYvectorFromFile(Const)
    %
    %   Input Arguments:
    %       Const
    %           A global struct containing program flow parameters and constants
    %       rhsfilename
    %           FEKO RHS vector filename (e.g. 'yagi.rhs')
    %
    %   Output Arguments:
    %       yVectors, a struct containing:
    %       values
    %           MxP set of Yrhs vectors for the MoM equation. Here,
    %           M is the number of (MoM) basis-functions read from the
    %           *.mat and *.str file
    %       numRhs
    %           The number of right hand sides (i.e. excitation vectors)
    %
    %   Description:
    %       Reads in a FEKO binary Yrhs vector file which contains either one
    %       or more excitations.
    %
    %   =======================
    %   Written by Danie Ludick on June 19, 2013.
    %   to account also for the fact that the *.rhs file can be generated
    %   with 32-bit or 64-bit FEKO
    %   Last updated on June 19, 2013.
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    %   Set constants related to which FEKO executables was used to generate the
    %   *.str file, i.e. either 32-bit and 64-bit FEKO executables. Please
    %   note that additional information on reading the *.mat files using
    %   FEKO can be found on the FEKO website:
    %   http://www.feko.info/support/helpcenter/how-to/how-to-read-the-.mat-.lud-.rhs-files-and-.str-files

    %   (Note, as of FEKO Suite 6.2 the BYTE_PADDING for the *.str file using both 32-bit and
    %   64-bit FEKO exectubles will be 4 Bytes) and the following settings
    %   should be fine for both 32-bit and 64-bit.
    PADDING_TYPE = 'int32'; % For 32-bit FEKO executables
    %PADDING_TYPE = 'int64' % For 64-bit FEKO executables (6.1.1 and earlier)

    error(nargchk(2,2,nargin));

    fid = fopen(rhsfilename,'r','l');

    if fid == -1
        error(['Error reading FEKO matrix file: ' rhsfilename]);
        mesage(Const,sprintf('Error read: %s',rhsfilename));
    end

    message_fc(Const,' ');
    message_fc(Const,'------------------------------------------------------------------------------------');
    message_fc(Const,sprintf('Reading *.rhs file from: %s',rhsfilename));

    MM = Const.numMoMbasis;
    yVectors.values = [];
    rawData = zeros(2*MM,1);
    precisionString = 'double'; % Yrhs is always double complex
    element_size = 8;
    numRhs = 1;
    while ~feof(fid)
        try
            % read the data for the column
            reclen = fread(fid,1,PADDING_TYPE);
            rawData = fread(fid,reclen/element_size,precisionString);
            reclen_end = fread(fid,1,PADDING_TYPE);
            % The start and end record length should match.
            if (reclen ~= reclen_end)
                message_fc(Const,'Invalid Record Length');
                error ('Invalid Record Length');
            end%if

            % If the record length is zero, then exit the loop
            if (reclen == 0)
              break;
            end%if
        catch
            % In this case an error has probably occurred.
            break;
        end
        if (~isempty(rawData))
            yVectors.values(1:MM,numRhs) = complex(rawData(1:2:2*MM),rawData(2:2:2*MM));
        end%if
        % Increment the number of RHS vectors
        numRhs = numRhs + 1;

	end%while
    fclose(fid);

    yVectors.numRhs = numRhs - 1;

    % Some debug information below
    %size(yVectors.values,1)
    %size(yVectors.values,2)

    message_fc(Const,sprintf('Read: %d excitations (number of Yrhs.)',yVectors.numRhs));
    message_fc(Const,sprintf('Finished processing the *.rhs file: %s',rhsfilename));