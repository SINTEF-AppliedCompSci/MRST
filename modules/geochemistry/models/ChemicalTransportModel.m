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
        fluidMat      % Matrix to compute, for a given component, the amount that is attached to the surface 
        surfMat       % Matrix to compute, for a given component, the amount that is contained in the fluid (water) 

    end


    methods

        function model = ChemicalTransportModel(G, rock, fluid, chemicalLogModel, varargin)

            model = model@WaterModel(G, rock, fluid, varargin{:});
            model.chemicalModel = chemicalLogModel;

            % Create a matrix of the components that are on surfaces
            chemmodel = model.chemicalModel;
            chemsys = chemmodel.chemicalSystem;
            nC        = chemsys.nC;
            nMC       = chemsys.nMC;
            CM        = chemsys.compositionMatrix;

            surfMaster  = logical(chemsys.surfMaster);
            surfComp    = sum(logical(CM(surfMaster, :)), 1);
            surfMult    = repmat(surfComp, nMC, 1);
            surfMatFlag = logical(CM.*surfMult);

            surfMat = zeros(nMC, nC);
            surfMat(surfMatFlag) = CM(surfMatFlag);
            model.surfMat = surfMat;
            
            fluidMat = zeros(nMC, nC);
            fluidMat(~surfMatFlag) = CM(~surfMatFlag);
            model.fluidMat = fluidMat;
            

        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)

            chemmodel = model.chemicalModel;
            
            [pVars, state] = prepStateForEquations(model, state);
               
            logSpecies  = chemmodel.getProp(state, 'logSpecies');
            logElements = chemmodel.getProp(state, 'logElements');
            
            species  = cellfun(@(x) exp(x), logSpecies, 'UniformOutput',false);
            elements = cellfun(@(x) exp(x), logElements, 'UniformOutput',false);

            state = chemmodel.setProp(state, 'species', species);
            state = chemmodel.setProp(state, 'elements', elements);
            
            [chem_eqs, chem_names, chem_types] = equationsChemicalLog(chemmodel, state);

            [tr_eqs, tr_names, tr_types] = equationsTransportComponents(model, ...
                                                              state0, state, ...
                                                              dt, ...
                                                              drivingForces);
            eqs = horzcat(tr_eqs, chem_eqs);
            names = {tr_names{:}, chem_names{:}};
            types = {tr_types{:}, chem_types{:}};

            problem = LinearizedProblem(eqs, types, names, pVars, state, dt);

        end

        function [variableNames, state] = prepStateForEquations(model, state)
            
            chemmodel = model.chemicalModel;
            chemsys = chemmodel.chemicalSystem;
            
            logSpeciesNames  = chemsys.logSpeciesNames;
            logElementNames  = chemsys.logElementNames;
            logSolidNames    = chemsys.logSolidNames;
            logGasNames      = chemsys.logGasNames;
            combinationNames = chemsys.combinationNames;
            logSurfActNames  = chemsys.logSurfaceActivityCoefficientNames;
            
            chemVariableNames = [ logSpeciesNames, logElementNames, logGasNames, logSolidNames, logSurfActNames, combinationNames];
            variableNames = {'pressure', chemVariableNames{:}};
            
            pressure = model.getProps(state, 'pressure');
            chemVariableValues = cell(1, numel(chemVariableNames));
            [chemVariableValues{:}] = chemmodel.getProps(state, chemVariableNames{:});
            variableValues = {pressure, chemVariableValues{:}};
            
            % initiate AD variables
            [variableValues{:}] = initVariablesADI(variableValues{:});
            chemVariableValues = {variableValues{2 : end}};
            
            % assign AD variables back to state.
            state = model.setProp(state, 'pressure', variableValues{1});
            state = chemmodel.setProps(state, chemVariableNames, chemVariableValues);            
            
        end
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments

            chemmodel = model.chemicalModel;
            chemsys   = chemmodel.chemicalSystem;
            
            vars = problem.primaryVariables;

            ind = false(size(vars));
            chemvars = {chemsys.logSpeciesNames{:}, chemsys.logElementNames{:},...
                        chemsys.logGasNames{:}, chemsys.logSolidNames{:},...
                        chemsys.logSurfaceActivityCoefficientNames{:},...
                        chemsys.combinationNames{:}}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(chemvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;

            chem_problem                  = problem;
            chem_problem.primaryVariables = vars(ind);
            chem_dx                       = dx(ind);
            
%             state = chemmodel.synclog(state);
            [state, chem_report] = chemmodel.updateState(state, chem_problem, ...
                                                         chem_dx, ...
                                                         drivingForces);


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
            
        end

        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
            % Note: Here component denotes in fact the fluid part of the master component.
            if isempty(force)
                return
            end

            chemmodel = model.chemicalModel;
            chemsys = chemmodel.chemicalSystem;
            
            ind = strcmpi(cname, chemsys.elementNames);
            if chemsys.surfMaster(ind)
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

        function [fn, index] = getVariableField(model, name, varargin)
                [fn, index] = getVariableField@WaterModel(model, name, varargin{:});
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
