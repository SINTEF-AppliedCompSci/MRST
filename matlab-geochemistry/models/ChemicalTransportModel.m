classdef ChemicalTransportModel < WaterModel
    % ChemicalTransportModel An object for the the simulation of chemical transport.
    %
    % SYNOPSIS:
    %  model = ChemicalTransportModel(G, rock, fluid, chemicalModel, varargin)
    %
    % DESCRIPTION:
    %   A class of WaterModel which can simulate advective transport of
    %   chemical through reactive media on arbitary domains. The solver
    %   uses a fully implicit numerical scheme with upwind discretization.
    %   The pressure, flow and chemistry system are fully coupled. 
    %
    % PARAMETERS:
    %   G - A geometry strucuture with the domain size and discretization 
    %
    %   rock - rock structure containing rock properties
    %
    %   fluid - fluid structure containing the properties of the fluid
    %           phase
    %
    %   chemicalModel - the chemical model object to be used during the
    %                   transport simulation
    %
    % RETURNS:
    %   model - the transport model object containing tools used internally for
    %           the solution of the coupled chemical, pressure and
    %           advection system. Useful fields include:
    %       model.fluidMat - matrix to determine the mass of each element
    %                       in the fluid phase
    %       model.plotIter - toggle the plotting of species and element
    %                        concentration in each newton step within a
    %                        timestep, [false]
    %       model.plotFinal - toggle the plotting of species and element
    %                         concentrations at the end of each time step, 
    %                         [false]
    %
    % EXAMPLES:
    %   %% define the domain
    %   G = cartGrid([100, 1, 1]);
    %   G = computeGeometry(G);
    % 
    %   %% Define the rock
    %   rock.perm = 1*darcy*ones(G.cells.num, 1);
    %   rock.poro = 0.5*ones(G.cells.num, 1);
    % 
    %   %% Define the fluid
    %   pRef = 0*barsa;
    %   fluid = initSimpleADIFluid('phases', 'W', 'mu', 1*centi*poise, 'rho', ...
    %                            1000*kilogram/meter^3, 'c', 1e-10, 'cR', 4e-10, ...
    %                            'pRef', pRef);
    % 
    %   %% Define the chemistry
    %   elements = {'O', 'H', 'Na*','Cl*'};
    % 
    %   species = {'H+*', 'OH-', 'Na+', 'H2O*', 'NaCl','Cl-'};
    % 
    %   reactions ={'H2O  = H+  + OH- ',      10^-14*mol/litre, ... 
    %             'NaCl = Na+ + Cl-',       10^1*mol/litre};
    % 
    %   % instantiate chemical model
    %   chemModel = ChemicalModel(elements, species, reactions);
    % 
    %   %% instantiate the transport model
    %   model = ChemicalTransportModel(G, rock, fluid, chemModel);
    %
    % SEE ALSO:
    %   'simulateScheduleAD', 'ChemicalModel'

    %{
    Copyright 2009-2017 SINTEF Digital, Mathematics & Cybernetics.

    This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

    MRST is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    MRST is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with MRST.  If not, see <http://www.gnu.org/licenses/>.
    %}
    
    properties
        
        chemicalModel % Chemical model for the chemistry
        chemical_fds % List of all the variable names for the chemistry
        transport_fds % List of all the variable names for the transport
        fluidMat % Matrix to compute, for a given component, the amount that is attached to the surface 
        surfMat % Matrix to compute, for a given component, the amount that is contained in the fluid (water) 
        plotIter % plot each iteration of the netwon solver
        plotFinal % plot final iteration of the newton solver
        currentTime % current time step in seconds

    end


    methods

        function model = ChemicalTransportModel(G, rock, fluid, chemicalLogModel, varargin)

            model = model@WaterModel(G, rock, fluid, varargin{:});
            model.chemicalModel = chemicalLogModel;
            model.chemical_fds = model.chemicalModel.getAllVarsNames();
            model.transport_fds = {'p', 'wellSol'};

            % Create a matrix of the components that are on surfaces
            chemModel = model.chemicalModel;
            nC        = chemModel.nC;
            nMC       = chemModel.nMC;
            CM        = chemModel.compositionMatrix;

            surfMaster  = logical(model.chemicalModel.surfMaster);
            surfComp    = sum(logical(CM(surfMaster, :)), 1);
            surfMult    = repmat(surfComp, nMC, 1);
            surfMatFlag = logical(CM.*surfMult);

            surfMat = zeros(nMC, nC);
            surfMat(surfMatFlag) = CM(surfMatFlag);
            model.surfMat = surfMat;
            
            fluidMat = zeros(nMC, nC);
            fluidMat(~surfMatFlag) = CM(~surfMatFlag);
            model.fluidMat = fluidMat;
            
            model.plotIter = false;
            model.plotFinal = false;


        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)

            [pVars, pressure, logComponents, logMasterComponents, combinations,...
                   logSaturationIndicies, logPartialPressures,...
                   logSurfaceAcitivityCoefficients] = prepStateForEquations(model, state);
               
            components = cellfun(@(x) exp(x), logComponents, 'UniformOutput',false);
            masterComponentss = cellfun(@(x) exp(x), logMasterComponents, 'UniformOutput',false);

            [chem_eqs, chem_names, chem_types] = equationsChemicalLog(model.chemicalModel, state, logComponents, logMasterComponents, combinations, ...
                                                       logPartialPressures, logSaturationIndicies,logSurfaceAcitivityCoefficients);


            [tr_eqs, tr_names, tr_types] = equationsTransportComponents(state0, ...
                                                              pressure, masterComponentss, ...
                                                              components,...
                                                              state, model, ...
                                                              dt, ...
                                                              drivingForces);
            eqs = horzcat(tr_eqs, chem_eqs );
            names = { tr_names{:},chem_names{:}};
            types = { tr_types{:},chem_types{:}};

            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end

        function [variableNames, pressure, logComponents, logMasterComponents, combinations,...
                   logSaturationIndicies, logPartialPressures,...
                   logSurfaceActivityCoefficients] = prepStateForEquations(model, state)
            
            chemModel = model.chemicalModel;

            logComponentNames       = chemModel.logSpeciesNames;
            logMasterComponentNames = chemModel.logElementNames;
            logSolidNames           = chemModel.logSolidNames;
            logGasNames             = chemModel.logGasNames;
            combinationNames        = chemModel.combinationNames;
            logSurfActNames         = chemModel.logSurfaceActivityCoefficientNames;
            
            
            variableNames = ['pressure', logComponentNames, logMasterComponentNames, logGasNames, logSolidNames, logSurfActNames, combinationNames];
            variableValues = cell(1, numel(variableNames));
            variableValues{1} = model.getProps(state, 'pressure');
            [variableValues{2:end}] = chemModel.getProps(state, variableNames{2:end});
            
            [variableValues{:}] = initVariablesADI(variableValues{:});
            
            logComponents        = cell(1, numel(logComponentNames));
            logMasterComponents  = cell(1, numel(logMasterComponentNames));
            logPartialPressures     = cell(1, numel(logGasNames));
            logSaturationIndicies   = cell(1, numel(logSolidNames));
            combinations   = cell(1, numel(combinationNames));
            logSurfaceActivityCoefficients = cell(1, numel(logSurfActNames));
            
            for i = 1 : numel(combinationNames)
                ind = strcmpi(combinationNames{i}, variableNames);
                combinations{i} = variableValues{ind};
            end

            for i = 1 : numel(logComponentNames)
                ind = strcmpi(logComponentNames{i}, variableNames);
                logComponents{i} = variableValues{ind};
            end

            for i = 1 : numel(logComponentNames)
                ind = strcmpi(logComponentNames{i}, variableNames);
                logComponents{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logMasterComponentNames)
                ind = strcmpi(logMasterComponentNames{i}, variableNames);
                logMasterComponents{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logGasNames)
                ind = strcmpi(logGasNames{i}, variableNames);
                logPartialPressures{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logSolidNames)
                ind = strcmpi(logSolidNames{i}, variableNames);
                logSaturationIndicies{i} = variableValues{ind};
            end
            
            for i = 1 : numel(logSurfActNames)
                ind = strcmpi(logSurfActNames{i}, variableNames);
                logSurfaceActivityCoefficients{i} = variableValues{ind};
            end
            
            
            ind = strcmpi('pressure', variableNames);
            pressure = variableValues{ind};
      
               
               
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments

            chemModel = model.chemicalModel;

            vars = problem.primaryVariables;

            ind = false(size(vars));
            chemvars = {chemModel.logSpeciesNames{:}, chemModel.logElementNames{:},...
                        chemModel.logGasNames{:}, chemModel.logSolidNames{:},...
                        chemModel.logSurfaceActivityCoefficientNames{:},...
                        chemModel.combinationNames{:}}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(chemvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;

            chem_problem                  = problem;
            chem_problem.primaryVariables = vars(ind);
            chem_dx                       = dx(ind);
            
%             state = chemModel.synclog(state);
            [state, chem_report] = chemModel.updateState(state, chem_problem, ...
                                                    chem_dx, drivingForces);


            ind = false(size(vars));
            fluidvars = {'pressure'}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;

            fluid_problem                  = problem;
            fluid_problem.primaryVariables = vars(ind);
            fluid_dx                       = dx(ind);

            [state, fluid_report] = updateState@WaterModel(model, state, ...
                                                           fluid_problem, ...
                                                           fluid_dx, ...
                                                           drivingForces);
            report = []; % no report for the moment.
            
            if model.plotIter
                h = findobj('tag', 'updateSpeciesfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updateSpeciesfig');
                    h = findobj('tag', 'updateSpeciesfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logSpecies);
                title('species - iteration');
                xlabel('cell index');
                ylabel('log_e species concentrations')
                legend(model.chemicalModel.speciesNames);
                drawnow;
                
                h = findobj('tag', 'updateMasterfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'updateMasterfig');
                    h = findobj('tag', 'updateMasterfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logElements);
                title('elements - iteration');
                xlabel('cell index');
                ylabel('log_e element concentrations')
                legend(model.chemicalModel.elementNames);
                drawnow;
                
                
            end

        end

        function [state, report] = updateAfterConvergence(model, state0, state, ...
                                                          dt, drivingForces) %#ok
            [state, report] = updateAfterConvergence@WaterModel(model, state0, ...
                                                              state, dt, drivingForces);
           
            
            if model.plotFinal 
                h = findobj('tag', 'convergedfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'convergedfig');
                    h = findobj('tag', 'convergedfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logSpecies);
                title('components - converged');
                xlabel('cell index');
                ylabel('log_e species concentrations')
                legend(model.chemicalModel.speciesNames);

                h = findobj('tag', 'convergedmasterfig');
                if isempty(h)
                    figure
                    set(gcf, 'tag', 'convergedmasterfig');
                    h = findobj('tag', 'convergedmasterfig');
                end
                set(0, 'currentfigure', h)
                clf
                plot(state.logElements);
                xlabel('cell index');
                ylabel('log_e element concentrations')
                title('master components - converged');
                legend(model.chemicalModel.elementNames);
                drawnow;
            end
            
        end


        function [eq, src] = addComponentContributions(model, cname, eq, ...
                                                       component, src, force)
            % Note: Here component denotes in fact the fluid part of the master component.
            if isempty(force)
                return
            end

            chemModel = model.chemicalModel;
            ind = strcmpi(cname, chemModel.elementNames);
            if chemModel.surfMaster(ind)
                return
            end

            c = model.getProp(force, cname);

            cells = src.sourceCells;
            qW = src.phaseMass{1}./model.fluid.rhoWS;

            isInj = qW > 0;

            qC = (isInj.*c + ~isInj.*component(cells)).*qW;

            eq(cells) = eq(cells) - qC;
            src.components{end+1} = qC;
        end

        function names = getComponentNames(model)
            names = model.chemicalModel.elementNames;
        end

        function [fn, index] = getVariableField(model, name)
            if ismember(name, model.chemical_fds)
                [fn, index] = model.chemicalModel.getVariableField(name);
            else
                [fn, index] = getVariableField@WaterModel(model, name);
            end
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@WaterModel(model);
            % Sources for chemical components
            forces.chemSrc = [];
            % chemSrc is a struc with fields
            % 'cells' : injection cells
            % 'comps' : compositions (for each cell)
            % 'rates' : injection rates (for each cell)
        end

    end
end
