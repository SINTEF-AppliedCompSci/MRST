classdef OilWaterFluidModel < TwoPhaseOilWaterModel
%
%
% SYNOPSIS:
%   model = OilWaterFluidModel(G, rock, fluid, varargin)
%
% DESCRIPTION: This model is derived from the fluid model TwoPhaseOilWaterModel
% and is added some few functionalities that are needed in the coupled solver.
%
% PARAMETERS:
%   G        - Grid structure
%   rock     - Rock structure
%   fluid    - Fluid structure
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO: MechFluidModel 
%


    properties
        % Primary variables that are used in the fluid system
        primaryVarNames;
    end

    methods
        function model = OilWaterFluidModel(G, rock, fluid, varargin)
        % Constructor
            model = model@TwoPhaseOilWaterModel(G, rock, fluid);
            model.primaryVarNames = {'pressure', 'sW'}; % well variables
                                                        % not included
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            error(['The OilWaterFluidModel is not meant for being called directly. ' ...
                   'It should be used coupled with a mechanical model, derived ' ...
                   'from MechFluidModel.'])
        end

        function varnames = getAllPrimaryVariables(model, state)
        % list all the primarty variables that are recognized and can be
        % handled by the model, used by the updateState member function.
            varnames = model.primaryVarNames;
            [~, wellVarNames, ~] = ...
                model.FacilityModel.getAllPrimaryVariables(state.wellSol);
            varnames = {varnames{:}, wellVarNames{:}};
        end

        function fds = getAllVarsNames(model)
        % list all the variables that are recognized and can be handled by the model
            fds = {'wellSol', 'pressure', 's', 'sW',  'sO', 'water', 'oil'};
        end

        function [state, report] = updateState(model, state, problem, dx, ...
                                               drivingForces)
            % updates the state variable that belongs to the model, that is the fluid
            % variables. The mechanical variables in the state will be updated
            % by the mechanical model, see MechanicMechModel.
            vars = problem.primaryVariables;
            ind = false(size(vars));
            fluidvars = model.getAllPrimaryVariables(state);

            [lia, loc] = ismember(fluidvars, vars);
            assert(all(lia), 'The primary variables are not set correctly.');
            ind(loc) = true;
            problem.primaryVariables = vars(ind);
            dx   = dx(ind);

            [state, report] = updateState@TwoPhaseOilWaterModel(model, state, ...
                                                              problem, dx, ...
                                                              drivingForces);

        end


    end

end
