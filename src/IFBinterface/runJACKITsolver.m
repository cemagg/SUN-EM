function [jack] = runJACKITsolver(Const, zMatrices, yVectors, xVectors, mbfs)
    %runJACKITsolver v0.2
    %   Date: 30.11.2013
    %   Usage:
    %       [jack] = runJACKITsolver(Const, zMatrices, yVectors, xVectors)
    %
    %   Input Arguments:
    %       Const
    %           A global struct, containing general data
    %       zMatrices
    %           The Z-matrices data
    %       yVectors
    %           The Yrhs-vector data
    %       xVectors
    %           The Xsol-vector data (i.e. MoM solution of FEKO)
    %       mbfs
    %           The MBFs (prim & sec.) generated in runMBFgenerator
    %
    %   Output Arguments:
    %       jack
    %           Struct containing the Jacobi Iterations solution and timing data
    %
    %   Description:
    %       Calculates the Jacobi Iterative Solver, as described in [1] and
    %       [2] and adapted for the DGFM in [3].
    %
    %   Assumptions:
    %        - All domains are the same size (i.e. contains the same number of unknowns)
    %
    %   TO-DO:
    %        - Currently only a single iteration is performed - we need to
    %          extend this to more than 1 iteration.
    %
    %   References:
    %   [1] Y. Brand, A. K. Skrivervik, J. R. Mosig, F. E. Gardoil, "New
    %       iterative integral equation technique for multilayered printed
    %       array antennas," in Mathematical Methods in Electromagnetic Theory,
    %       Kharkov, Ukraine, Jun, 1998, pp. 615-617.
    %   [2] Brand, Yan, Anja K. Skrivervik, and Juan R. Mosig. "An iterative scheme solution
    %       for the analysis of printed arrays." Microwave and Optical Technology Letters
    %       16.2 (1997): 106-115.
    %   [3] A. C. Polycarpou, "Evaluation of stationary block iterative
    %       techniques for the solution of finite arrays using fe-bi method
    %       and domain decomposition," in Proc. European Conference on
    %       Antennas and Propagation (EuCAP), Nice, France, Nov. 2006,
    %       pp. 1-6
    %   [4] D.J. Ludick, R. Maaskant, et. al, "Efficient Analysis of Large
    %       Irregular Antenna Arrays using the Domain Green's Function Method",
    %       Special Issue on Antennas Propagation Submission, 2013 (pending)
    %   =======================
    %   Written by Danie Ludick on June 24, 2013.
    %   Last updated on November 30, 2013.
    %   EMSS-SA (Pty) Ltd
    %   Email: dludick@emss.co.za

    error(nargchk(5,5,nargin));

    % Formating based on where the routine is called from
    if (Const.runJACKITfromDGFM)
        message_fc(Const,sprintf('  Running Jacobi Iterative solver'));
        if (Const.no_mutual_coupling_array)
            message_fc(Const,sprintf('  (*** Mutual coupling between array elements ignored ***)'));
        end%if
    else
        message_fc(Const,' ');
        message_fc(Const,'------------------------------------------------------------------------------------');
        message_fc(Const,sprintf('Running Jacobi Iterative solver'));
        if (Const.no_mutual_coupling_array)
            message_fc(Const,sprintf('(*** Mutual coupling between array elements ignored ***)'));
        end%if
    end%if

    if (Const.no_mutual_coupling_array)
        % If now mutual coupling should be accounted in the array environment, then just switch off calculating secondary MBFs
        Const.calcSecMBFs = false;
    end%if

    % Start timing
    tic

    % Initialise the return values
    jack  = [];
    jack.name = 'jack';
    Nmom        = Const.numMoMbasis;            % Total number of basis functions for whole problem
    Nngf        = Const.numNGFbasis;            % Number of basis functions for NGF domain
    NtotArray   = Const.numArraybasis;          % Total number of basis functions for the array
    Ndom        = Const.numMoMbasisPerElement;  % Number of basis functions per array element
    numArrayEls = Const.numArrayElements;       % The number of array elements
    numSols     = xVectors.numSols;             % The number of reference solutions
    jack.numSols = numSols;                     % Calculate a solution for each configuration
    jack.Isol   = complex(zeros(NtotArray,numSols));   % The Jacobi method is used to calculate the array solution
    % Extract the solutions from which to start and end our CBFM process
    solStart = Const.solStart;
    solEnd   = Const.solEnd;

    if (Const.JACKITiterScheme >= 0)
        jack.numIter = Const.JACKITiterScheme;
        jack.relError = zeros(jack.numIter+1,numSols); % +1 because we start the Jacobi iterations at n=0 below
    else
        jack.numIter = 1;
        jack.relError = zeros(numSols);
    end%if

    % TO-DO: See FEKDDM-4.3: Check convergence of the Jacobi Iterative Solver
    if (checkJACKITconvergence(Const, zMatrices) ~= true)
       error (['No convergence possible for Jacobi Iter. Solver']);
       message_fc(Const,sprintf('No convergence possible for Jacobi Iter. Solver'));
    end

    % TO-DO: See issue FEKDDM-10: We have not yet added support for multiple solutions here for
    % Const.JACKITiterScheme values other than -1 and also comment in the driver routine.

    if (Const.JACKITiterScheme >= 0)

        % TO-DO: This currently does not support the use of multiple
        % solution configurations.
        if (numSols > 1)
            error (['Cannot calculate Jacobi Iterations with Scheme > 0 for multiple solution configurations']);
            message_fc(Const,sprintf('Cannot calculate Jacobi Iterations with Scheme > 0 for multiple solutions'));
        end%if

        % -- Multiple iterations (See FEKDDM-4.2), with the number of iterations specified
        % Implemented here, Eq. 18 [3] - J = Sum_{n=0}^{K}{-(Zon)^(-1)Zoff}^n*J0
        % Note, J0 is the primary MBFs that was calculated earlier and K is
        % limited to jack.numIter.

        % Add some progress output here
        h = waitbar(0,'Starting Jacobi Iterations ...');

        ZonInv = complex(zeros(Ndom,Ndom));
        Zoff = complex(zeros(Ndom,Ndom));
        % Because we are working with the same domains, we can extract and
        % calculate only the inverse of one diagonal element. This can then
        % be used subsequently.
        ZonInv = inv(zMatrices.values((1-1)*Ndom+1:1*Ndom,(1-1)*Ndom+1:1*Ndom));
        for n=0:jack.numIter
            % Extract first the self-interaction / diagonal matrix here (Zon)
            % See issue FEKDDM-4.2: Not required for each domain, as the
            % domains are identical.
            % Zon = zMatrices.values((p-1)*Ndom+1:m*Ndom,(p-1)*Ndom+1:p*Ndom);
            for p=1:numArrayEls
                if (n==0)
                    % Take into account the addition of the primary MBF for domain p (i.e. the 0th order current component)
                    for prim_ii=1:mbfs.numPrimMBFs(p)
                        % See issue FEKDDM-4.4: Check whether this is indeed the correct manner to treat multiple
                        % primary MBFs induced on a domain, e.g. for multiple ports.
                        jack.Isol((p-1)*Ndom+1:p*Ndom) = jack.Isol((p-1)*Ndom+1:p*Ndom) + mbfs.PrimIsol(:,prim_ii,p);
                    end%for
                    % Skip to next loop iteration
                    continue;
                end%if
                for m=1:numArrayEls
                    if (m==p)
                        % Self-term already accounted for.
                        continue;
                    end%if
                    % Extract the coupling matrix, Zoff
                    Zoff = zMatrices.values((p-1)*Ndom+1:p*Ndom,(m-1)*Ndom+1:m*Ndom);
                    for prim_ii=1:mbfs.numPrimMBFs(m)
                        % See issue FEKDDM-4.4: Check whether this is indeed the correct manner to treat multiple primary MBFs induced
                        % on a domain, e.g. for multiple ports.
                        jack.Isol((p-1)*Ndom+1:p*Ndom) = jack.Isol((p-1)*Ndom+1:p*Ndom) + (-ZonInv*Zoff)^n * mbfs.PrimIsol(:,prim_ii,m);
                    end%for
                end
            end%for p=1:numArrayEls
            waitbar(n/jack.numIter,h,sprintf('Calculating Jacobi Iter. solution: %.2f%% complete',100*(n/jack.numIter)))
            % Compare the solution obtained here for this iteration, with that obtained by FEKO that was stored in xVectors.values
            jack.relError(n+1) = calculateErrorNormPercentage(xVectors.values, jack.Isol);
        end%for n=1:jack.numIter

        close(h); % Close the progress bar
        plotData1(Const, [0:n],jack.relError, 'Iteration (N)','Rel. Error Percentage', 'The Relative Err. Percentage for Jacobi Iteration compared to the FEKO MoM solution', 'jacobiDipoleArrayTightCoupling')

    elseif (Const.JACKITiterScheme == -1)

        % See issue FEKDDM-10: We added now support for multiple solution
        % configurations. Each of the solution configurations has given rise to
        % a suitable weightVectors (based on the excitation-law, exact current,
        % etc. for that specific solution). We need to repeat the DGFM active
        % impedance matrix calculation of each array element for each solution
        % configuration.
        for solNum = solStart:solEnd
            % -- Use the primary and secondary MBFs generated in the generic MBF routine
            for p=1:numArrayEls
                % Add the contribution of the primary MBFs
                for prim_ii=1:mbfs.numPrimMBFs(p,solNum)
                    % Note, assumed here for mbfs.PrimIsol, is that each of the domains only have 1
                    % primary MBF. The last term in the following assignment has to be changed to:
                    % mbfs.PrimIsol(:,m,p)
                    jack.Isol((p-1)*Ndom+1:p*Ndom,solNum) = jack.Isol((p-1)*Ndom+1:p*Ndom,solNum) + mbfs.PrimIsol(:,prim_ii,p,solNum);
                end%for n=1:numArrayEls
                if (~Const.calcSecMBFs)
                    % If there are no Secondary MBFs, then continue with the
                    % next array element (only primary MBFs then taken into
                    % account)
                    continue;
                end%if

                %  Add the contribution of the secondary MBFs
                secIndx = 0;
                for sec_ii=1:mbfs.numSecMBFs(p,solNum)
                    % See FEKDDM-4.1 and also comment in runMBFgenerator: The negative sign has already been
                    % accounted for when the secondary MBFs are generated.
                    jack.Isol((p-1)*Ndom+1:p*Ndom,solNum) = jack.Isol((p-1)*Ndom+1:p*Ndom,solNum) + mbfs.SecIsol(:,sec_ii,p,solNum);
                end%for
            end%for

            % See issue FEKDDM-6.2: If the NGF-enhanced DGFM is used, then we
            % compare only the current that is on the dynamic domain, i.e. the
            % finite array. This routine might also be called from the IFB-DGFM solver (therefore used now Const.domAoffset instead of Nngf)
            jack.relError(solNum) = calculateErrorNormPercentage(xVectors.values(Const.domAoffset+1:Nmom,solNum), jack.Isol(:,solNum));
        end%for solNum = 1:numSols

    elseif (Const.JACKITiterScheme == -2)
        % -- Multiple iterations (See FEKDDM-4.2), iterate until convergence is reached
        % (Just use a WHILE-loop here, or we can merge this with the above)
        error(['This iteration scheme is not yet supported for Jacobi Iterative Solver']);
        mesage(Const,sprintf('This iteration scheme is not yet supported for Jacobi Iterative Solver'));
    end%if

    % End timing
    jack.solTime = toc;

    % Write the JACKIT solution to a ASCII str file, so that it can be read
    % again by FEKO (for plotting in POSTFEKO) - only if requested
    if (isequal(and([0 0 0 0 1 0 0],Const.FEKDDMwriteFEKOstrfile),[0 0 0 0 1 0 0]))
        writeSolToFile(Const, jack);
    end%if

    % Do output for each solution configuration
    for solNum = solStart:solEnd

        % Formating based on where the routine is called from
        if (Const.runJACKITfromDGFM)
            % Note, use here jack.numIter+1 because we include also the 0th iteration
            message_fc(Const,sprintf('  Finished Jacobi Iterative solver in %f sec. with %d iterations for Sol. %d of %d', ...
                jack.solTime, jack.numIter+1,solNum,numSols));
            message_fc(Const,sprintf('  Rel. error norm. compared to FEKO MoM sol. %f percent for Sol. %d of %d',...
                jack.relError(solNum),solNum,numSols));
        else
            message_fc(Const,sprintf('Finished Jacobi Iterative solver in %f sec. with %d iterations for Sol. %d of %d',...
                jack.solTime, jack.numIter+1,solNum,numSols));
            message_fc(Const,sprintf('Rel. error norm. compared to FEKO MoM sol. %f percent for Sol. %d of %d',...
                jack.relError(solNum),solNum,numSols));
        end%if
    end