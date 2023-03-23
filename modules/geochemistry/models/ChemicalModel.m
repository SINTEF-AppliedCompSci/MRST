classdef ChemicalModel < PhysicalModel

    properties
        
        chemicalSystem
        
        % Solver parameters
        plotIterations            % toggle plotting of iteration information true or false
        nonLinearMinIterations    % minimum number of iterations for the nonlinear solver
        nonLinearMaxIterations    % maximum number of iterations for the nonlinear solver  
        nonLinearTolerance        % tolerance of the residual of the nonlinear system
        linearTolerance           % tolerance of the residual of the linear system, for backslash
        linearMaxIterations       % maximum number of iterations for the linear solver
        
    end


    methods

        function model = ChemicalModel(chemsys, varargin)

            model = model@PhysicalModel([], varargin{:});
            model.chemicalSystem = chemsys;

            model.plotIterations            = false;
            model.nonLinearMinIterations    = 1; 
            model.nonLinearMaxIterations    = 25;  
            model.nonLinearTolerance        = 1e-12;
            model.linearMaxIterations       = 25;
            model.linearTolerance           = 1e-8;
            
        end


        %%
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)

            chemsys = model.chemicalSystem;
            [logcomps, logmasterComps] = prepStateForEquations(model, state);

            [eqs, names, types] = equationsChemicalLog(logcomps, logmasterComps, ...
                                                       model);

            primaryVariables = chemsys.logSpeciesNames;
            problem = LinearizedProblem(eqs, types, names, primaryVariables, state, dt);

        end


        %%
        function [fn, index] = getVariableField(model, name, varargin)

            chemsys = model.chemicalSystem;
            varfound = false;
         

            while ~varfound


                if strcmpi(name, 'partialPressures')
                    varfound = true;
                    fn = 'partialPressures';
                    index = ':';
                    break
                end

                ind = strcmpi(name, chemsys.partialPressureNames);
                if any(ind)
                    varfound = true;
                    fn = 'partialPressures';
                    index = find(ind);
                    break
                end

                if strcmpi(name, 'fluidVolumeFraction')
                    varfound = true;
                    fn = 'fluidVolumeFraction';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'logFluidVolumeFraction')
                    varfound = true;
                    fn = 'logFluidVolumeFraction';
                    index = ':';
                    break
                end

                 if strcmpi(name, 'matrixVolumeFraction')
                    varfound = true;
                    fn = 'matrixVolumeFraction';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'logMatrixVolumeFraction')
                    varfound = true;
                    fn = 'logMatrixVolumeFraction';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'species')
                    varfound = true;
                    fn = 'species';
                    index = ':';
                    break
                end

                if strcmpi(name, 'chargeBalance')
                    varfound = true;
                    fn = 'chargeBalance';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'elements')
                    varfound = true;
                    fn = 'elements';
                    index = ':';
                    break
                end

                if strcmpi(name, 'logSpecies')
                    varfound = true;
                    fn = 'logSpecies';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'activities')
                    varfound = true;
                    fn = 'activities';
                    index = ':';
                    break
                end

                if strcmpi(name, 'logElements')
                    varfound = true;
                    fn = 'logElements';
                    index = ':';
                    break
                end

                if strcmpi(name, 'combinationComponents')
                    varfound = true;
                    fn = 'combinationComponents';
                    index = ':';
                    break
                end
                 
                ind = strcmpi(name, chemsys.combinationNames);
                if any(ind)
                    varfound = true;
                    fn = 'combinationComponents';
                    index = find(ind);
                    break
                end
              
                
                ind = strcmpi(name, chemsys.speciesNames);
                if any(ind)
                    varfound = true;
                    fn = 'species';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'aqueousConcentrations')
                    varfound = true;
                    fn = 'aqueousConcentrations';
                    index = ':';
                    break
                end

                ind = strcmpi(name, chemsys.aqueousConcentrationNames);
                if any(ind)
                    varfound = true;
                    fn = 'aqueousConcentrations';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'surfaceConcentrations')
                    varfound = true;
                    fn = 'surfaceConcentrations';
                    index = ':';
                    break
                end

                ind = strcmpi(name, chemsys.surfaceConcentrationNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceConcentrations';
                    index = find(ind);
                    break
                end

                ind = strcmpi(name, chemsys.gasNames);
                if any(ind)
                    varfound = true;
                    fn = 'partialPressures';
                    index = find(ind);
                    break
                end

                if strcmpi(name, 'saturationIndicies')
                    varfound = true;
                    fn = 'saturationIndicies';
                    index = ':';
                    break
                end

                ind = strcmpi(name, chemsys.solidNames);
                if any(ind)
                    varfound = true;
                    fn = 'saturationIndicies';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, chemsys.solidDensityNames);
                if any(ind)
                    varfound = true;
                    fn = 'solidDensities';
                    index = find(ind);
                    break
                end
                
                 if strcmpi(name, 'solidDensities')
                    varfound = true;
                    fn = 'solidDensities';
                    index = ':';
                    break
                end
                
                if strcmpi(name, 'logPartialPressures')
                    varfound = true;
                    fn = 'logPartialPressures';
                    index = ':';
                    break
                end

                ind = strcmpi(name, chemsys.logGasNames);
                if any(ind)
                    varfound = true;
                    fn = 'logPartialPressures';
                    index = find(ind);
                    break
                end

                if strcmpi(name, 'logSaturationIndicies')
                    varfound = true;
                    fn = 'logSaturationIndicies';
                    index = ':';
                    break
                end

                ind = strcmpi(name, chemsys.logSolidNames);
                if any(ind)
                    varfound = true;
                    fn = 'logSaturationIndicies';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, chemsys.activityNames);
                if any(ind)
                    varfound = true;
                    fn = 'activities';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, chemsys.logSpeciesNames);
                if any(ind)
                    varfound = true;
                    fn = 'logSpecies';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, chemsys.elementNames);
                if any(ind)
                    varfound = true;
                    fn = 'elements';
                    index = find(ind);
                    break
                end

                ind = strcmpi(name, chemsys.logElementNames);
                if any(ind)
                    varfound = true;
                    fn = 'logElements';
                    index = find(ind);
                    break
                end

                ind = strcmpi(name, 'temperature');
                if any(ind)
                    varfound = true;
                    fn = 'temperature';
                    index = find(ind);
                    break
                end
              
                
                ind = strcmpi(name, chemsys.surfaceChargeNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceCharges';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, chemsys.surfacePotentialNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfacePotentials';
                    index = find(ind);
                    break
                end
                
                ind = strcmpi(name, chemsys.surfaceActivityCoefficientNames);
                if any(ind)
                    varfound = true;
                    fn = 'surfaceActivityCoefficients';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'surfaceActivityCoefficients')
                    varfound = true;
                    fn = 'surfaceActivityCoefficients';
                    index = ':';
                    break
                end
                
                ind = strcmpi(name, chemsys.logSurfaceActivityCoefficientNames);
                if any(ind)
                    varfound = true;
                    fn = 'logSurfaceActivityCoefficients';
                    index = find(ind);
                    break
                end
                
                if strcmpi(name, 'logSurfaceActivityCoefficients')
                    varfound = true;
                    fn = 'logSurfaceActivityCoefficients';
                    index = ':';
                    break
                end
                
                [fn, index] = getVariableField@PhysicalModel(model, name, varargin{:});
                varfound = true;
            end

        end

        %%
        function fds = getAllVarsNames(model)
            chemsys = model.chemicalSystem;
            fds = {'species', 'elements', 'logSpecies', ...
                   'logElements', chemsys.speciesNames{:}, ...
                   chemsys.elementNames{:}, chemsys.logSpeciesNames{:}, ...
                   chemsys.logElementNames{:}, 'logFluidVolumeFraction', chemsys.logSurfaceActivityCoefficientNames{:}, chemsys.logGasNames{:}, chemsys.logSolidNames{:}};
        end

        %%
        function [state, report, model] = initState(model, userInput, varargin)
        % initState solves the chemical system using the specified input
        %
        % SYNOPSIS:
        %  [state, report, model] = initState(model, userInput, varargin)
        %
        % DESCRIPTION:
        %   A function used to solve the chemical system generated by
        %   ChemicalModel according to constraints given in userInput. 
        %
        % REQUIRED PARAMETERS:
        %   userInput - A matrix containing the values of elements
        %               or species specified as inputs during the generation of the
        %               chemical system using ChemicalModel. userInput must have as
        %               many columns as there are specified knowns, which is also equal
        %               to the number of entries in elements. 
        %           
        %
        % KEYWORD ARGUMENTS:
        %   'temperature' - 'key'/value pair representing the
        %                   temperature of the chemical system. Can be a scalar, or a
        %                   vector with the same length as userInput. If temperature is
        %                   not provided 298 K is used. 
        %          
        %    	state = model.initState(userInput, 'temp', 305)
        %
        %   'chargeBalance' - Toggles the enforcment of strict
        %                     charge balance. A master component must be specified as a 
        %                     'key'/value pair. If charge balance is not satisfied with the
        %                     given constraint the solver will adjust the total concentration
        %                     of the specified element to enforce charge balance. The 
        %                     specified element must be one of the elements provided as an
        %                     input. By defualt charge balance is not strictly enforced. 
        %                     This make the chemical solver more robust. 
        %
        %    	state = model.initState(userInput, 'chargeBalance', 'Cl')
        %
        %
        % RETURNS:
        %   state - A structure containing the solution to the
        %           chemical system. This includes state.elements and
        %           state.species. state.elements contains the total
        %           element concentrations and surface functional group
        %           concentrations. state.species contains the species
        %           concentrations. Values can be retreived from state using the
        %           model.getProp() function.
        %
        %   	OH = chem.getProp(state, 'OH-');
        %
        %
        %   report - The nonlinear solver report for the chemical system. 
        %
        %   model - The updated chemical model object.
        %
        % SEE ALSO:
        %   'ChemicalModel', 'getProps'
        

            chemsys = model.chemicalSystem;
            % parse inputs to initState
            p = inputParser;
            
            valFun = @(x) ischar(x);
            p.addParameter('chargeBalance', 'nochargebalance', valFun);
            p.addParameter('temperature', 298, @isnumeric);

            
            p.parse(varargin{:})
            
            if strcmpi('nochargebalance',p.Results.chargeBalance)
               
            elseif any(ismember({'e','e-'},p.Results.chargeBalance))
                error('Electron concentration, e or e-, can not be used for charge balance.')
                
            elseif sum(ismember(chemsys.elementNames,p.Results.chargeBalance))==0
                warning('Using any quantity other than a total element concentration as a charge balance can yield unexpected results.')
                
            end
            
                
            nI = size(userInput,1);
            
            % grab temperature
            nTemp = size(p.Results.temperature, 1);
            if nTemp == 1
                state.temperature = repmat(p.Results.temperature, nI, 1);
            else
                assert(nTemp == nI, 'The number of cells in state.temperature do not correspond to the size of userInput. Repmat the input vector if you want to do a temperature sweep.');
                state.temperature = p.Results.temperature(:);
            end

            
            givenTest = any(strcmpi(p.Results.chargeBalance, horzcat(chemsys.inputNames,'nochargebalance')));
            if ~givenTest
                warning 'Using a total element concentration whos values are not given (marked with "*") for charge balance may result in unexpected results.';
            end
            
            chargeBalance = ~strcmpi(p.Results.chargeBalance,'nochargebalance');

            chemicalInputModel       = ChemicalInputModel(chemsys);
            compositionReactionModel = CompositionReactionModel(chemsys);
            
            inputNames = chemicalInputModel.inputNames;

            inSize = size(userInput);
            if ~isempty(chemsys.surfInfo);
                k = numel(chemsys.surfInfo.master);
            else
                k = 0;
            end
            
            assert(inSize(2) == chemsys.nMC-k, ['For the specified chemical system the input constraint must have ' ...
                                num2str(chemsys.nMC-k) ' columns.']);
            
            ncells = size(userInput,1);
            state.elements = mol/litre*ones(ncells, chemsys.nMC);
            state.species  = mol/litre*ones(ncells, chemsys.nC);
            
            if ~isfield(state, 'combinationComponents') && chemsys.nLC > 0
                state.combinationComponents = ones(ncells, chemsys.nLC);
            end 
            
            if ~isfield(state, 'partialPressures') && chemsys.nG > 0
                state.partialPressures = 0.1*ones(ncells, chemsys.nG);
            end 
            
            if ~isfield(state, 'saturationIndicies') && chemsys.nS > 0
                state.saturationIndicies = 0.1*ones(ncells, chemsys.nS);
            end 

            if ~isfield(state, 'surfaceActivityCoefficients') && chemsys.nP > 0
                state.surfaceActivityCoefficients = 0.1*ones(ncells, chemsys.nP);
            end 

            state = model.syncLog(state);
            state = model.validateState(state);

            % put the value of the surface functional groups into the input
            % vector
            call = 0;
            for i = 1 : chemsys.nMC
                
                if regexpi(inputNames{i}, '>')
                    sInd = strcmpi(inputNames{i}, chemsys.surfInfo.master);
                    d = chemsys.surfInfo.d{sInd};
                    s = chemsys.surfInfo.s{sInd};
                    a = chemsys.surfInfo.a{sInd};
                    in{i} = d.*s.*a;
                    call = call + 1;
                else
                    in{i} = userInput(:,i-call);
                end
                
            end
            
            % make sure input is the correct size
            givenSize = cellfun(@(x) size(x,1), in);
            mSize = max(givenSize);
            assert(any(givenSize == mSize | givenSize == 1), 'Input size must be consistent between input and surface parameters');
           
            % set their values inside of state
            for i = 1 : chemsys.nMC
            	state = model.setProp(state, inputNames{i}, in{i}.*ones(mSize,1));
            end
            
            state = model.syncLog(state);
            
            % create initial guess
            fprintf('Computing initial guess...\n')
            [state, ~, report] = compositionReactionModel.solveChemicalState(state);
            
            % solve chemical system
                
            if chargeBalance
                fprintf('Solving the chemical system with strict charge balance...\n');
                state0 = state;
                chargeBalanceModel = ChargeBalanceModel(chemsys);
                chargeBalanceModel.CVC = p.Results.chargeBalance;
                [state, ~, report] = chargeBalanceModel.solveChemicalState(state);
                
                if ~report.Converged
                    state = state0;
                    warning('Charge balance not acheived with given inputs, consider choosing a different component for charge compensation. Attempting solution without charge balance.');
                end
            end
            if ~report.Converged || ~chargeBalance
                fprintf('Solving chemical system...\n')
                [state, ~, report] = chemicalInputModel.solveChemicalState(state);
                if ~report.Converged
                    warning('Solver did not converge, use results with caution.');
                end
            end
        end

        %%
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments
            state0 = state;
            [state, report] = updateState@PhysicalModel(model, state, problem, ...
                                                        dx, drivingForces);

            state = model.syncFromLog(state);
            state = updateChemicalModel(model, problem, state, state0 );
            state = model.syncLog(state);
            

            if model.plotIterations
                h = findobj('tag', 'updatechemfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updatechemfig');
                    h = findobj('tag', 'updatechemfig');
                end
                set(0, 'currentfigure', h)
                clf
                if size(state.species, 1) ~= 1
                    plot(log10(state.species*litre/mol));
                    title('species (chemistry step)');
                    legend(model.speciesNames);
                end

                h = findobj('tag', 'chemmastercomp');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'chemmastercomp');
                    h = findobj('tag', 'chemmastercomp');
                end
                set(0, 'currentfigure', h)
                clf
                if size(state.species, 1) ~= 1
                    plot(log10(state.elements*litre/mol));
                    title('master species (chemistry step)');
                    legend(model.elementNames);
                end
                drawnow;
            end

        end

        %%
        function state = validateState(model, state)
            state = validateState@PhysicalModel(model, state);

            chemsys = model.chemicalSystem;
            
            prop = 'species';
            fn = model.getVariableField(prop);
            assert(isfield(state, fn), ['property ', prop, ' must be supplied']);
            
            prop = 'elements';
            fn = model.getVariableField(prop);
            assert(isfield(state, fn), ['property ', prop, ' must be supplied']);
            
            prop = 'combinationComponents';
            if chemsys.nLC > 0
                fn = model.getVariableField(prop);
                assert(isfield(state, fn), ['property ', prop, ' must be supplied']);
            end 
            
            prop = 'partialPressures';
            if chemsys.nG > 0
                fn = model.getVariableField(prop);
                assert(isfield(state, fn), ['property ', prop, ' must be supplied']);
            end 
            
            prop = 'saturationIndicies';
            if chemsys.nS > 0
                fn = model.getVariableField(prop);
                assert(isfield(state, fn), ['property ', prop, ' must be supplied']);
            end 
            
            prop = 'surfaceActivityCoefficients';
            if chemsys.nP > 0
                fn = model.getVariableField(prop);
                assert(isfield(state, fn), ['property ', prop, ' must be supplied']);
            end 

            state = model.syncLog(state);
        end

        %%
        function state = syncLog(model, state)
            state.logElements = log(state.elements);
            state.logSpecies  = log(state.species);
            
            if isfield(state, 'saturationIndicies');
                state.logSaturationIndicies = log(state.saturationIndicies);
            end
            if isfield(state, 'partialPressures');
                state.logPartialPressures = log(state.partialPressures);
            end
            if isfield(state, 'surfaceActivityCoefficients');
                state.logSurfaceActivityCoefficients = log(state.surfaceActivityCoefficients);
            end
        end

        %%
        function state = syncFromLog(model, state)
            state.elements = exp(state.logElements);
            state.species       = exp(state.logSpecies);
            
            if isfield(state, 'logSaturationIndicies');
                state.saturationIndicies     = exp(state.logSaturationIndicies);
            end
            if isfield(state, 'logPartialPressures');
                state.partialPressures       = exp(state.logPartialPressures);
            end
            if isfield(state, 'logSurfaceActivityCoefficients');
                state.surfaceActivityCoefficients       = exp(state.logSurfaceActivityCoefficients);
            end
        end


        %%
        function [logComps, logMasterComps] = prepStateForEquations(model, state)

            chemsys = model.chemicalSystem;
            logComps = cell(chemsys.nC, 1);
            [logComps{:}] = model.getProps(state, chemsys.logSpeciesNames{:});

            logMasterComps = cell(chemsys.nMC, 1);
            [logMasterComps{:}] = model.getProps(state, chemsys.logElementNames{:});

            [logComps{:}] = initVariablesADI(logComps{:});

        end
        
        function p = getPropAsCell(model, state, fd)
            p = model.getProp(state, fd);
            if ~iscell(p)
                nrow = size(p, 2);
                p = arrayfun(@(i) p(:, i), 1 : nrow, 'uniformoutput', false);
            end
        end
        

        %%
        function [state, failure, report] = solveChemicalState(model, inputstate)
        % inputstate contains the input and the initial guess.

            inputstate = model.validateState(inputstate); % in particular,
                                                          % updates the log
                                                          % variables if necessary.

            solver = NonLinearSolver('maxIterations', model.maxIterations);
            dt = 0; % dummy timestep
            drivingForces = []; % drivingForces;
            inputstate0 = inputstate;

            [state, failure, report] = solveMinistep(solver, model, inputstate, ...
                                                     inputstate0, dt, ...
                                                     drivingForces);

        end

        %%
        function [state, model] = updateActivities(model, state)
        % updateAcitivities updates the acitivity of each aqueous species
        % using the extended Davies equaiton and adds the values to the 
        % field state.activities.
        %
        % SYNOPSIS:
        %  [state, chem] = updateActivities(chem, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState. 
        %          
        %
        % RETURNS:
        %   state - The state structure containing the additional field 
        %           state.acitivities which holds the acitivity of each aqueous species 
        %           in units of mol/meter^3. Activities can be retrieved using the 
        %           getProp command by calling for the species name prepended by 'a.'
        %
        %   chem - updated chemical object including the field chem.activityNames 
        %          for use in calling the activities using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.updateActivities(state);
        %   aH2O = chem.getProp(state, 'aH2O');
        %   pH = -log10(chem.getProp, 'aH+');
        %
        % SEE ALSO:
        %   'updateChargeBalance'
        
            assert( nargout == 2, ['Output argument to updateActivities must include the state variable and the chemical model object']);
            [state, model] = computeActivity(model, state);
            
        end
        
        %%
        function [state, model] = updateChargeBalance(model, state)
        % updateChargeBalance updates the residual of the aqueous charge
        % balance equation and adds the values to the field state.chargeBalance.
        %
        % SYNOPSIS:
        %  [state, chem] = updateChargeBalance(chem, state)
        %
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState.
        %
        % RETURNS:
        %   state - The state structure containing additional field state.chargeBalance
        %           which holds the aqueous charge balance residual in units of fraction of total
        %           charge. The charge balance can be retrieved using the getProps
        %           command by calling for the variables 'chargeBalance'.
        %
        %   chem  - updated chemical object including for use in calling the
        %           charge balance residual using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.updateChargeBalance(state);
        %   charge = chem.getProp(state, 'chargeBalance');
        %
        % SEE ALSO:
        %   'updateActivities'        

            assert( nargout == 2, ['Output argument to updateChargeBalance must include the state variable and the chemical model object']);
            [state, model] = computeChargeBalance(model, state);
            
        end
        
        %%
        function [state, model] = updateSurfacePotentials(model, state)
        % updateSurfacePotentials updates the potential of each layer of each
        % surface and adds the values to the field state.surfacePotentials.
        %
        % SYNOPSIS:
        %  [state, chem] = updateSurfacePotentials(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by model.initState.
        %          
        % RETURNS:
        %   state - The state structure containing the additional field state.surfacePotentials of
        %           which holds the potential of each layer of
        %           each surface in Volts. Values can be retrieved using
        %           the getProps command by calling for the surface functional
        %           group name, followed by '_Psi_' followed by the layer
        %           number (0, 1, 2) or for the constant capacitance model just
        %           use '_Psi'.
        %
        %   chem - updated chemical object including the field chem.surfacePotentialNames 
        %           for use in calling the surface potentials using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.updateSurfacePotentials(state);
        %   [pot1, pot2] = chem.getProps(state, '>SiO_Psi_0', '>SiO_Psi_1');
        %
        % SEE ALSO:
        %   'updateSurfaceCharges'
        
        
            assert( nargout == 2, ['Output argument to updateSurfacePotentials must include the state variable and the chemical model object']);
            [state, model] = computeSurfacePotential(model, state);
            
        end
        
        %%
        function [state, model] = updateSurfaceCharges(model, state)
        % updateSurfaceCharges updates the charge of each layer of each
        % surface and adds the values to the field state.surfaceCharges.
        %
        % SYNOPSIS:
        %  [state, chem] = updateSurfaceCharges(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState.
        %          
        % RETURNS:
        %   state - The state structure with the additional field 
        %           state.surfaceCharges containing the charge of each layer of
        %           each surface in C/meter^2. Values can be retrieved using
        %           the getProps command by calling for the surface functional
        %           group name, followed by '_sig_' followed by the layer
        %           number (0, 1, 2) or for the constant capacitance model just
        %           use '_sig'.
        %
        %   chem - updated chemical object including the field chem.surfaceChargeNames 
        %           for use in calling the surface charges using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.updateSurfaceCharges(state);
        %   charge0 = chem.getProp(state, '>SiO_sig_0');
        %
        % SEE ALSO:
        %   'updateSurfacePotentials'
        

            assert( nargout == 2, ['Output argument to updateSurfaceCharges must include the state variable and the chemical model object']);                
            [state, model] = computeSurfaceCharge(model, state);
        end
        
        %%
        function [state, model] = updateAqueousConcentrations(model, state)
        % updateAqueousConcentrations updates the total concentration of each 
        % element in the aqueous phase. This excludes, gas, solid and
        % surface bound species. 
        %
        % SYNOPSIS:
        %  [state, chem] = updateAqueousConcentrations(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by chem.initState.
        %          
        % OUTPUTS:
        %   state - The state structure containing the additional field 
        %           state.aqueousConcentrations which hold the total aqueous
        %           concentration of each element in the aqueous phase. The values can be retrieved using the
        %           getProps function by asking for the element name followed
        %           by '(aq)'.
        %
        %   chem - updated chemical object including the field chem.surfaceChargeNames 
        %          for use in calling the surface charges using getProp/s.
        %
        % EXAMPLE:
        %   [state, chem] = chem.updateAqueousConcentrations(state);
        %   Na = chem.getProp(state, 'Na(aq)');
        %   H = chem.getProp(state, 'H(aq)');
        %
        % SEE ALSO:
        %   'updateSurfaceConcentrations'
        
        
            assert( nargout == 2, ['Output argument to updateAqueousConcentrations must include the state variable and the chemical model object']);
            [state, model] = computeAqueousConcentrations(model, state);
            
        end
        
        %%
        function [state, model] = updateSurfaceConcentrations(model, state)
        % updateSurfaceConcentrations updates the total concentration of each 
        % element bound to surfaces.
        %
        % SYNOPSIS:
        %  [state, chem] = updateAqueousConcentrations(chem, state)
        %
        % REQUIRED PARAMETERS:
        %   state - the state variable produced by model.initState.
        %          
        % OUTPUTS:
        %   state - The state structure containing the additional field 
        %           state.surfaceConcentrations which holds the surface
        %           concentration of each element. The values can be retrieved using the
        %           getProps function by asking for the element name followed
        %           by '_surf'.
        %
        %   chem - updated chemical object now includes the field 
        %           chem.surfaceConcentractionNames for use in  calling the 
        %           surface concetrations using getProp/s. 
        %
        % EXAMPLE:
        %   [state, chem] = chem.updateSurfaceConcentrations(state);
        %   [Na, H] = chem.getProps(state, 'Na(surf)', 'H(surf)');
        %
        % SEE ALSO:
        %   'updateSurfaceConcentrations'
        

            assert( nargout == 2, ['Output argument to updateSurfaceConcentrations must include the state variable and the chemical model object']);
            [state, model] = computeSurfaceConcentrations(model, state);
            
        end
    end

end
