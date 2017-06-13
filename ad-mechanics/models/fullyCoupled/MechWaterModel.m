classdef MechWaterModel < MechFluidModel
%
%
% SYNOPSIS:
%   model = MechWaterModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION: Model for coupled mechanical fluid simulation. The fluid model
% is a single phase model.
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
%{
Copyright 2009-2017 SINTEF ICT, Applied Mathematics.

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

    methods
        function model = MechWaterModel(G, rock, fluid, mech_problem, varargin)

            model = model@MechFluidModel(G, rock, fluid, mech_problem, ...
                                         varargin{:});

        end

        function fluidModel = setupFluidModel(model)
            fluidModel = WaterFluidModel(model.G, model.rock, ...
                                               model.fluid);
        end


        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, wellSol, xd] = model.getProps(state, 'pressure', 'wellSol', ...
                                                     'xd');
            % Properties at previous timestep
            [p0, xd0] = model.getProps(state0, 'pressure', 'xd');

            %Initialization of primary variables

            fluidModel = model.fluidModel; % shortcuts
            mechModel  = model.mechModel; % shortcuts


            [wellVars, wellVarNames, wellMap] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            [p, wellVars{:}, xd] = initVariablesADI(p, wellVars{:}, xd);

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd);

            [w_eqs, w_eqsnames, w_eqstypes, state] = equationsWaterMech(state0, ...
                                                              p, wellVars, ...
                                                              state, fluidModel, ...
                                                              dt, mechTerm, ...
                                                              drivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            eqs = horzcat(w_eqs, mech_eqs);
            names = {w_eqsnames{:}, mech_eqsnames{:}};
            types = {w_eqstypes{:}, mech_eqstypes{:}};

            primaryVars = {'pressure', wellVarNames{:}, 'xd'};
            % make sure that the primary variables defined here match with
            % those of WaterFluidModel and MechanicalModel.

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function [mechTerm, fluidp] = computeCouplingTerms(model, p0, ...
                                                              xd0, p, xd)

            opmech = model.mechModel.operators;
            fluidp = p;
            mechTerm.new = opmech.div*xd;
            mechTerm.old = opmech.div*xd0;

        end


    end
end
