classdef ChemicalTransportModel < WaterModel

    properties
        % Chemical model for the chemistry
        chemicalModel
        % List of all the variable names for the chemistry
        chemical_fds;
        % List of all the variable names for the transport
        transport_fds;
    end


    methods

        function model = ChemicalTransportModel(G, rock, fluid, chemicalModel, varargin)

            model = model@WaterModel(G, rock, fluid, varargin{:});
            model.chemicalModel = chemicalModel;
            model.chemical_fds = model.chemicalModel.getAllVarsNames();
            model.transport_fds = {'p', 'wellSol'};
            
            % TODO check if there is no collisions in the name.

        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            [p, wellSol] = model.getProps(state, 'pressure', 'wellSol');
            
            
            chemModel = model.chemicalModel;
            % chemical equations

            compNames = chemModel.CompNames;
            masterCompNames = chemModel.MasterCompNames;

            comps = cell(1, numel(compNames));
            mastercomps = cell(1, numel(masterCompNames));
            
            [comps{:}] = model.getProps(state, compNames{:});
            [mastercomps{:}] = model.getProps(state, masterCompNames{:});
            
            [p, comps{:}, mastercomps{:}] = initVariablesADI(p, comps{:}, ...
                                                             mastercomps{:});
            
            [chem_eqs, chem_names, chem_types] = equationsChemical(comps, mastercomps, ...
                                                    chemModel);

            [tr_eqs, tr_names, tr_types] = equationsTransportComponents(state0, ...
                                                              p, mastercomps, ...
                                                              comps,...
                                                              state, model, ...
                                                              dt, ...
                                                              drivingForces);

            primaryVars = {'pressure', chemModel.CompNames{:}, chemModel.MasterCompNames{:}};
            eqs = horzcat(tr_eqs, chem_eqs );
            names = { tr_names{:},chem_names{:}};
            types = { tr_types{:},chem_types{:}};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces) %#ok
        % Update state based on Newton increments

            chemModel = model.chemicalModel;
            
            vars = problem.primaryVariables;

            ind = false(size(vars));
            chemvars = {chemModel.CompNames{:}, chemModel.MasterCompNames{:}}; % the chemical primary variables, see getEquations
            [lia, loc] = ismember(chemvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;
            
            chem_problem                  = problem;
            chem_problem.primaryVariables = vars(ind);
            chem_dx                       = dx(ind);

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
        end

        function [eq, src] = addComponentContributions(model, cname, eq, component, src, force)
        % For a given component conservation equation, compute and add in
        % source terms for a specific source/bc where the fluxes have
        % already been computed.
        %
        % INPUT:
        %
        % model  - (Base class, automatic)
        %
        % cname  - Name of the component. Must be a property known to the
        %          model itself through getProp/getVariableField.
        %
        % eq     - Equation where the source terms are to be added. Should
        %          be one value per cell in the simulation grid (model.G)
        %          so that the src.sourceCells is meaningful.
        %
        % component - Cell-wise values of the component in question. Used
        %          for outflow source terms only.
        %
        % src    - Source struct containing fields for fluxes etc. Should
        %          be constructed from force and the current reservoir
        %          state by computeSourcesAndBoundaryConditionsAD.
        %
        % force  - Force struct used to produce src. Should contain the
        %          field defining the component in question, so that the
        %          inflow of the component through the boundary condition
        %          or source terms can accurately by estimated.
            if isempty(force)
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
            names = model.chemicalModel.MasterCompNames;
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
