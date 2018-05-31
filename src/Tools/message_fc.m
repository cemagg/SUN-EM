function message_fc(Const,Text)

% 10-June-2013 Danie Ludick, EM System & Software (Pty) Ltd (FEKO), dludick@emss.co.za
% 2017.12.12: We need to make sure we restore again the current working directory!

% 2018-01-02: Renamed now the routine to message_fc.m, as message.m results in problems
% when using MATLAB R2017b.

cwd = pwd();

cd(Const.ProjectPath);
cd(Const.OutputDirName);
fid=fopen('SimOutput.log','a');
if fid==-1
    disp('Error: cannot open or create file SimOutput.log')
else
    C=fix(clock);
    fprintf(fid,['%d-%d-%d, %d:%d:%d => ' Text '\n'],C(3),C(2),C(1),C(4),C(5),C(6));
    disp(Text);
    fclose(fid);
end

% Flush the output buffer. Note, the manner in which this is done, is
% different in octave and MATLAB
if (Const.is_octave)
    fflush(stdout);
else
    drawnow('update');
end%if

% 2017.12.12: We need to make sure we restore again the current working directory!
cd(cwd);

return
