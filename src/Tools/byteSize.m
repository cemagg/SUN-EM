function [memStr] = byteSize(var)
    %byteSize v1.0
    %   Date: 21.06.2013
    %   Usage:
    %       memStr = byteSize(x)
    %
    %   Input Arguments:
    %       Const
    %           The variable whos memory usage is to be determined (note,
    %           for structs, use the getfield property to extract the
    %           struct member)
    %
    %  Output Arguments:
    %       memStr
    %           String that contains the memory usage (in b, kb, Mb, etc.)
    %
    %   Description:
    %       Returns the memory usage of the provide variable
    %
    %   =======================
    %   Written by Danie Ludick on June 21, 2013
    %   Last updated on June 21, 2013.
    %   EM Systems & Software (Pty) Ltd.
    %   Email: dludick.emss.co.za

    error(nargchk(1,1,nargin));
    
	s = whos('var');
	memStr = Bytes2str(s.bytes);

function str = Bytes2str(NumBytes)
% BYTES2STR Private function to take integer bytes and convert it to
% scale-appropriate size.

	scale = floor(log(NumBytes)/log(1024));
	switch scale
        case 0
            str = [sprintf('%.0f',NumBytes) ' b'];
        case 1
            str = [sprintf('%.2f',NumBytes/(1024)) ' kb'];
        case 2
            str = [sprintf('%.2f',NumBytes/(1024^2)) ' Mb'];
        case 3
            str = [sprintf('%.2f',NumBytes/(1024^3)) ' Gb'];
        case 4
            str = [sprintf('%.2f',NumBytes/(1024^4)) ' Tb'];
        case -inf
            % Size occasionally returned as zero (eg some Java objects).
            str = 'Not Available';
        otherwise
           str = 'Over a petabyte!!!';
	end%switch
