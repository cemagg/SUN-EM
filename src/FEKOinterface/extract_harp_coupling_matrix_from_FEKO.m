function [Zpq] = extract_harp_coupling_matrix_from_FEKO(Const, Solver_setup, r_T)
    %extract_harp_coupling_matrix_from_FEKO v1.0
    %   Date: 2017.12.12
    %   Usage:
    %       [Zpq] = extract_harp_coupling_matrix_from_FEKO(Const, Solver_setup, r_T)
    %
    %   Input Arguments:
    %       Const: A global struct containing:
    %         harp_mbf_array_FEKOprefilename
    %           FEKO *.pre filename (e.g. 'strip_dipole_2_el_array.pre') that is the
    %           *.pre file for a finite array that will be used to extract a 
    %           coupling matrix (used in the HARP MBF generation of secondaries)
    %           This file essentially contains 2 x array elements. The first is positioned
    %           at the origin, and the second at an arbitrary position specified by r_T.
    %       Solver_setup:
    %           Details about the solver setup, e.g. number of basis
    %           functions, etc.
    %       r_T : The position of the secondary element
    %   Output Arguments:
    %       Zpq:
    %           The coupling submatrix calculated using FEKO (based on the array elements and
    %           and positions as specified in the *.pre and *.xml files respectively)
    %
    %   Description:
    %       Reads in a FEKO *.pre and *.xml file and extracts the coupling submatrix from the *.mat file
    %
    %   =======================
    %   Written by Danie Ludick on 2017.12.12
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za
    
    error(nargchk(3,3,nargin));

    LOCAL_DEBUG = false; % Some local debug info

    % Save the current working directory - we are going to change into the project directory to process the 
    % file
    cwd = pwd();
    cd(Const.ProjectPath);

    % A debug variable (used to count the secondary MBF #, so that we can
    % append this to the XML file
    persistent secondary_mbf_index; % Initialise to zero
    
    % Upon the first call we need to initialize the variable.
    if isempty(secondary_mbf_index)
        secondary_mbf_index = 0;
    end
    
    % Open the *.pre file - we do this, just to check that it is there.
    fid_pre_file = fopen(Const.harp_mbf_array_FEKOprefilename,'r');
    if fid_pre_file == -1
        message_fc(Const,sprintf('Error reading HARP MBF array FEKO *.pre file: %s',Const.harp_mbf_array_FEKOprefilename));
        error(['Error reading FEKO *.pre file: ' Const.harp_mbf_array_FEKOprefilename]);
    end
    fclose(fid_pre_file); % Close it again - no need to modify anything.

    % Open the *.XML file
    fid_xml_file = fopen(Const.harp_mbf_array_FEKOxmlfilename,'w+');
    if fid_xml_file == -1
        message_fc(Const,sprintf('Error reading HARP MBF array FEKO *.xml file: %s',Const.harp_mbf_array_FEKOxmlfilename));
        error(['Error reading FEKO *.xml file: ' Const.harp_mbf_array_FEKOxmlfilename]);
    end    

    message_fc(Const,sprintf('      Generating coupling matrix for HARP using the FEKO files:'));
    message_fc(Const,sprintf('        *.pre file: %s',Const.harp_mbf_array_FEKOprefilename) )
    message_fc(Const,sprintf('        *.xml file: %s',Const.harp_mbf_array_FEKOxmlfilename) )

    % Update the array position in the XML file
    message_fc(Const,sprintf('        Updating second array position to: (x,y,z) = (%f,%f,%f)',r_T(1),r_T(2),r_T(3)));


    % ----------------------------------------------------
    % Update the XML file for the array positions
    % ----------------------------------------------------
 
    % First write out the header information
    fprintf(fid_xml_file,'<?xml version="1.0" encoding="utf-8"?>\n');
    fprintf(fid_xml_file,'<ArrayDistributionMatrix xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xmlns:xsd="http://www.w3.org/2001/XMLSchema" Version="2.0.0.0" DesignFrequency="%d" FrequencyUnit="Hz" CoordinateUnit="Meters" xmlns="http://www.antennamagus.com/schemas/ArrayDistributionMatrix.xsd"> <Elements>\n',Const.harps_freq_MHz*1E6);
    fprintf(fid_xml_file,'            <Element Name="%d" X="%f" Y="%f" Z="%f" Magnitude="1.0" Phase="0.0" PhiRotation="0" ThetaRotation="0" GammaRotation="0" />\n',1,0,0,0);
    fprintf(fid_xml_file,'            <Element Name="%d" X="%f" Y="%f" Z="%f" Magnitude="1.0" Phase="0.0" PhiRotation="0" ThetaRotation="0" GammaRotation="0" />\n',2,r_T(1),r_T(2),r_T(3));
    fprintf(fid_xml_file,'    </Elements>\n');
    fprintf(fid_xml_file,'</ArrayDistributionMatrix>');

    % Close again the XML file
    fclose(fid_xml_file);

    % ----------------------------------------------------
    % Using the updated XML file, run feko_mkl (assume it is in the $PATH)
    % ----------------------------------------------------
    if (ispc || isunix)
        % ----------
        % Run PREFEKO
        % ----------
        message_fc(Const,sprintf('        Running PREFEKO'))
        % Call prefeko for the pre file (that now contains the updated array element position)
        [stat, sysout] = system(sprintf('prefeko %s', Const.harp_mbf_array_FEKOprefilename));
        %sysout
        if (stat ~= 0)
            sysout % Write the FEKO output to the screen
            message_fc(Const,sprintf('Error running prefeko for: %s', Const.harp_mbf_array_FEKOprefilename));
            message_fc(Const, num2str(stat));
            error(['Error running prefeko for: ' Const.harp_mbf_array_FEKOprefilename]);
        end%if

        % ----------
        % Run FEKO (sequential for now)
        % ----------
        message_fc(Const,sprintf('        Running FEKO_MKL (Sequential)'))

        % Run feko_mkl directly (which will also beforehand run *.prefeko) - 
        % remove the *.pre extension before doing so
        fekfilename = regexprep(Const.harp_mbf_array_FEKOprefilename,'.pre','');

        % Call feko_mkl for the fek file (that now contains the updated array element position)
        [stat, sysout] = system(sprintf('feko_mkl %s', fekfilename));
        if (stat ~= 0)
            sysout % Write the FEKO output to the screen
            message_fc(Const,sprintf('Error running feko_mkl for: %s', fekfilename));
            message_fc(Const, num2str(stat));
            error(['Error running feko_mkl for: ' fekfilename]);            
        end%if

    else
        error(['Error: Running feko_mkl for HARP solver only possible on Windows or Linux']);
        mesage_fc(Const, sprintf('Error: Running feko_mkl for HARP solver only possible on Windows Linux'));
    end%if

    % --------------------
    % Extract now the *.mat file from the results
    % --------------------

    % Set a flag to indicate we are reading it for this solution step of HARP
    % Flag no longer needed, as we pass the *.mat filename to the readFEKOZMatrixFromFile routine
    %Const.harp_coupling_matrix_calcActive = true;
    Const.harp_mbf_array_FEKOmatfilename = regexprep(Const.harp_mbf_array_FEKOprefilename,'.pre','.mat');
    [zMatrices_harp] = readFEKOZMatrixFromFile(Const, Const.harp_mbf_array_FEKOmatfilename);
    %Const.harp_coupling_matrix_calcActive = false;

    % Extract here the coupling matrix
    num_bfs_per_el = Solver_setup.mom_basis_functions_per_array_element; 
    ObservRWGs = [1:num_bfs_per_el];
    SourceRWGs = [num_bfs_per_el+1:2*num_bfs_per_el];

    % NOTE: We assume that we are only working with a single freq here
    Zpq = calcZmn(Const, zMatrices_harp, 1, 1, 1, ObservRWGs, SourceRWGs);    
    
    if (LOCAL_DEBUG)
        % Increment our secondary MBF index used below in the debug output
        secondary_mbf_index = secondary_mbf_index + 1;
        
        % For debugging purposes, it is good to save each of the XML files that
        % are generated for the MBFs, to ensure that we are generating them
        % correctly        
        old_filename = Const.harp_mbf_array_FEKOxmlfilename;
        new_filename = sprintf('Sec_MBF_%d_%s',secondary_mbf_index,Const.harp_mbf_array_FEKOxmlfilename);        
        [err, msg] = movefile(old_filename, new_filename);
    end
    
    % Before exiting, go back to the working directory that was active before we changed into the
    % project directory for this routine:
    cd(cwd);
        
    % Exit the routine now
    return