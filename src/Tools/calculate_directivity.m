function calculate_directivity(Const, Solver_setup, Solution, theta_list_degrees, phi_list_degrees)
    %calculate_directivity
    %   Date: 2018.05.23
    %   Usage:
    %       calculate_directivity(Const, Solver_setup, Solution, theta_list, phi_list)
    %
    %   Input Arguments:
    %       Const: A global struct containing:
    %       Solver_setup (struct):
    %           Solver setup data like frequency range, rwg indices, triangle data, etc.
    %       Solution (struct):
    %           Struct containing the expansion coefficients, as well as other solver data
    %           e.g. timing info
    %       theta_list_degrees, phi_list_degrees:
    %           Vector of theta and phi co-ordinates (in degrees)
    %
    %   Output Arguments:
    %       None
    %
    %   Description:
    %       Calculates the directivity patterns at the various angular co-ordinates
    %
    %   =======================
    %   Written by Danie Ludick on 2018.05.18
    %   Stellenbosch University
    %   Email: dludick.sun.ac.za

