classdef MechBlackOilModel < MechFluidModel
%
%
% SYNOPSIS:
%   model = MechBlackOilModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION:
%   Model for fully coupled mechanical fluid simulation. The fluid
%   model is a blackoil model.
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
%   MechOilWaterModel, MechWaterModel, MechFluidModel

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
        function model = MechBlackOilModel(G, rock, fluid, mech_problem, varargin)

            model = model@MechFluidModel(G, rock, fluid, mech_problem, ...
                                         varargin{:});

        end

        function fluidModel = setupFluidModel(model)
            fluidModel = BlackOilFluidModel(model.G, model.rock, ...
                                                 model.fluid);
            fluidModel.disgas = true;
        end


        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only
            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, sW, sG, rs, rv, wellSol, xd] = model.getProps(state, 'pressure', ...
                                                                     'water', ...
                                                                     'gas', ...
                                                                     'rs', ...
                                                                     'rv', ...
                                                                     'wellSol', ...
                                                                     'xd');
            % Properties at previous timestep
            [p0, sW0, sG0, rs0, rv0, xd0] = model.getProps(state0, 'pressure', ...
                                                                   'water', ...
                                                                   'gas', 'rs', ...
                                                                   'rv', ...
                                                                   'xd');

            %Initialization of primary variables

            fluidModel = model.fluidModel; % shortcuts
            mechModel  = model.mechModel; % shortcuts

            st  = fluidModel.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
            st0 = fluidModel.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);
            if fluidModel.disgas || fluidModel.vapoil
                % X is either Rs, Rv or Sg, depending on each cell's saturation status
                x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
                gvar = 'x';
            else
                x = sG;
                gvar = 'sG';
            end

            [wellVars, wellVarNames, wellMap] = ...
                fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);

            if ~opt.resOnly,
                % define primary varible x and initialize
                [p, sW, x, wellVars{:}, xd] = initVariablesADI(p, sW, x, ...
                                                               wellVars{:}, ...
                                                               xd);
            end

            [mechTerm, fluidp] = computeCouplingTerms(model, p0, xd0, p, xd, ...
                                                      state0.u, state.u);

            [bo_eqs, bo_eqsnames, bo_eqstypes, state] = equationsBlackOilMech(state0, ...
                                                              st0, p, sW, x, ...
                                                              rs, rv, st, ...
                                                              wellVars, state, ...
                                                              fluidModel, dt, ...
                                                              mechTerm, ...
                                                              drivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsPoroMechanics(xd, ...
                                                              mechModel, ...
                                                              fluidp);

            eqs = horzcat(bo_eqs, mech_eqs);
            names = {bo_eqsnames{:}, mech_eqsnames{:}};
            types = {bo_eqstypes{:}, mech_eqstypes{:}};

            primaryVars = {'pressure', 'sW', gvar, wellVarNames{:}, 'xd'};
            % make sure that the primary variables defined here match with
            % those of BlackOilFluidModel and MechanicalModel.

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
           % why we divide by the cell's volumes.
           
        end
    end
end
