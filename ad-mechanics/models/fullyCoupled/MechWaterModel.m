classdef MechWaterModel < MechFluidModel
%
%
% SYNOPSIS:
%   model = MechWaterModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION:
%   Model for fully coupled mechanical fluid simulation. The fluid model is
%   a single phase model.
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
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   MechBlackOilModel, MechOilWaterModel, MechFluidModel

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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
                                         'water', true, varargin{:});
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

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd, ...
                                                      state0.u, state.u);

            [w_eqs, w_eqsnames, w_eqstypes, state] = equationsWaterMech(p0, state0, ...
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

        function [problem, state] = getAdjointEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', true, ...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                     'wellSol' , ...
                                                     'xd');
            % Properties at previous timestep
            [p0, wellSol0, xd0] = model.getProps(state0, 'pressure', ...
                                                         'wellSol' , ...
                                                         'xd');

            fluidModel = model.fluidModel; % shortcuts
            mechModel  = model.mechModel;  % shortcuts

            [wellVars, wellVarNames, wellMap] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            [wellVars0, ~, ~] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol0);

            if opt.reverseMode
                [p0, wellVars0{:}, xd0] = initVariablesADI(p0, wellVars0{:}, xd0);
            else
                [p, wellVars{:}, xd] = initVariablesADI(p, wellVars{:}, xd);
            end

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd, ...
                                                      state0.u, state.u);

            [ow_eqs, ow_eqsnames, ow_eqstypes, state] = equationsWaterMech( ...
                                                              p0,  state0, ...
                                                              p,  wellVars, state, ...
                                                              fluidModel, ...
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

            if opt.reverseMode
                primaryVars = {'pressure', 'xd'};
            else
                primaryVars = {'pressure', 'xd', wellVarNames{:}};
            end

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end
        
        
        function [mechTerm, fluidp] = ...
               computeCouplingTerms(model, p0, xd0, p, xd, u0, u)

           G = model.G;
           op = model.mechModel.operators;
           fluidp = p;
           
           if isa(xd, 'ADI')
              u = double2ADI(u, xd);
              u(~op.isdirdofs) = xd;  % to ensure correct derivatives
           end

           mechTerm.new = (op.ovol_div*u)./(G.cells.volumes);
           mechTerm.old = (op.ovol_div*u0)./(G.cells.volumes);
            
           % Note that the opmech.div returns the divergence integrated over cells. That is
           % why we divide by the cell's volumes

        end

    end
end
