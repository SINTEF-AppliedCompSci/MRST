classdef MechOilWaterModel < MechFluidModel
%
%
% SYNOPSIS:
%   model = MechOilWaterModel(G, rock, fluid, mech_problem, ...
%
% DESCRIPTION: Model for coupled mechanical fluid simulation. The fluid model
% is a two phase oil water model.
%
% PARAMETERS:
%   G            - grid structure
%   rock         - rock structure
%   fluid        - fluid structure
%   mech_problem - Structure that contains the mechanical parameters of the system
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%
% SEE ALSO:
%


    methods
        function model = MechOilWaterModel(G, rock, fluid, mech_problem, ...
                                           varargin)

            model = model@MechFluidModel(G, rock, fluid, mech_problem, ...
                                         varargin{:});

        end

        function fluidModel = setupFluidModel(model)
            fluidModel = OilWaterFluidModel(model.G, model.rock, model.fluid);
        end


        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, sW, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                         'water', ...
                                                         'wellSol', ...
                                                         'xd');
            % Properties at previous timestep
            [p0, sW0, xd0] = model.getProps(state0, 'pressure', ...
                                                    'water', ...
                                                    'xd');

            %Initialization of primary variables

            fluidModel = model.fluidModel; % shortcuts
            mechModel  = model.mechModel; % shortcuts

            [wellVars, wellVarNames, wellMap] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly,
                % define primary varible x and initialize
                [p, sW, wellVars{:}, xd] = initVariablesADI(p, sW, wellVars{:}, ...
                                                            xd);
            end

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd);

            [ow_eqs, ow_eqsnames, ow_eqstypes, state] = equationsOilWaterMech(state0, ...
                                                              p, sW, wellVars, ...
                                                              state, fluidModel, ...
                                                              dt, mechTerm, ...
                                                              drivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            eqs = horzcat(ow_eqs, mech_eqs);
            names = {ow_eqsnames{:}, mech_eqsnames{:}};
            types = {ow_eqstypes{:}, mech_eqstypes{:}};

            primaryVars = {'pressure', 'sW',  wellVarNames{:}, 'xd'};
            % make sure that the primary variables defined here match with
            % those of OilWaterFluidModel and MechanicalModel.

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, ...
                                                           xd0, p, xd)

            G = model.G;
            op = model.mechModel.operators;
            fluidp = p;
            mechTerm.new = (op.div*xd)./(G.cells.volumes);;
            mechTerm.old = (op.div*xd)./(G.cells.volumes);
            % Note that the opmech.div returns the divergence integrated over cells. That is
            % why we divide by the cell's volumes

        end


    end
end
