function centrePoint = lineCentre(p1, p2)
    %lineCentre
    %   Usage:
    %       centrePoint = lineCentre(p1, p2)
    %
    %   Input Arguments:
    %       p1, p2
    %           The two points in cartesian co-ordinates forming the end-nodes of the line
    %   Output Arguments:
    %       centrePoint
    %           The cartesian co-ordinate of the centre-points between p1 and p2
    %
    %   Description:
    %       Determines the centre-point of a line, given the
    %              end-points in cartesian coordinates
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.04
    %   Stellenbosch University
    %   Email: dludick@sun.ac.za

    centrePoint = (p1 + p2) .* 0.5;
