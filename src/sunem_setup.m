function [Const] = sunem_setup(Const)
    %sunem_setup
    %   Usage:
    %           [Const] = sunem_setup(Const, 'exampleA')
    %
    %   Input Arguments:
    %       Const
    %           A global program flow and configuration struct
    %   Output Arguments:
    %       Const
    %           Updated global struct (Const) containing the correct file
    %           paths, etc.
    %
    %   Description:
    %       Sets up the environment and also sets the correct version number
    %
    %   =======================
    %   Written by Danie Ludick on 2015-04-10
    %   Email: dludick@sun.ac.za

    % --------------------------------------------------------------------------------------------------
    % Set the version Major.Minor.Patch (see CHANGELOG.md for details)
    % --------------------------------------------------------------------------------------------------
    Const.MajorVersion = 1;  % Major version
    Const.MinorVersion = 35;  % Minor version

    Const.Year = 2019;
    Const.Month = 08;
    Const.Day = 12;

    % --------------------------------------------------------------------------------------------------
    % Set the Project path directories, and also add the tools, and interfaces to the Path
    % --------------------------------------------------------------------------------------------------
    cd(Const.ProjectPath);
    cd ..; cd ..;
    MainPath=pwd;
    cd(Const.ProjectPath);
    warning off MATLAB:MKDIR:DirectoryExists;
    mkdir(Const.OutputDirName);
    mkdir([Const.OutputDirName '/Temp']);

    % === General Tools  ===
    addpath([MainPath '/src/Tools']);

    % === FEKO Interface ===
    addpath([MainPath '/src/FEKOinterface']);

    % === FEKO Interface (binaries) ===
    addpath([MainPath '/bin/']);

    % === MoM Solver Interface ===
    addpath([MainPath '/src/MoMinterface']);

    % === MoM Solver Interface as implemented in [DBD2011] ===
    addpath([MainPath '/src/dBMoMinterface']);
    
    % === MBF Solver Interface ===
    addpath([MainPath '/src/MBFinterface']);
    
    % === DGFM Solver Interface ===
    addpath([MainPath '/src/DGFMinterface']);

    % === Iterative Solver Interface ===
    addpath([MainPath '/src/IFBinterface']);

    % === CMA Interface ===
    addpath([MainPath '/src/CMAinterface']);
    
    % ==  Set the executables path % ==
    if (ispc)
        Const.ExecPath = sprintf('%s\\bin\\',MainPath);
    else
        % Assume Unix like environment (e.g. Linux or Mac)
        Const.ExecPath = sprintf('%s//bin',MainPath);
    end%if

    % Set the project name
    [pathstr, Const.ProjectName, ext] = fileparts(pwd);

    % Check whether we are working in Octave or in MATLAB and save this as a constant
    result = exist('octave_config_info');
    if (result)
        Const.is_octave = true;
    else
        Const.is_octave = false;
    end%if

    % --------------------------------------------------------------------------------------------------
    % Output to screen and log-file
    % --------------------------------------------------------------------------------------------------
    message_fc(Const,' ');
    message_fc(Const,'==================================================================================');
    message_fc(Const,'======                          SUN-EM                                      ======');
    message_fc(Const,sprintf('======                 Version %d.%d Date: %d-%d-%d                         ======',...
        Const.MajorVersion,Const.MinorVersion, Const.Year,Const.Month,Const.Day));
    % TO-DO: Danie, add also here supported *.fek file formats
    message_fc(Const,'======                                                                      ======');
    message_fc(Const,'======      CEM toolkit and solvers from Stellenbosch University            ======');
    message_fc(Const,'======                      Copyright 2018                                  ======');
    if (Const.debug)
        message_fc(Const,'======             **********    DEBUG RUN   **********                     ======');
    end%if
    message_fc(Const,'==================================================================================');
    message_fc(Const,' ');

    % --------------------------------------------------------------------------------------------------
    % Some additional settings
    % --------------------------------------------------------------------------------------------------
    format long         % Set the format to long output

