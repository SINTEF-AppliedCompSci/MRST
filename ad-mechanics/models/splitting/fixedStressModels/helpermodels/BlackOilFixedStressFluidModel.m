classdef BlackOilFixedStressFluidModel < ThreePhaseBlackOilModel
%
%
% SYNOPSIS:
%   model = BlackOilFixedStressFluidModel(G, rock, fluid, varargin)
%
% DESCRIPTION:
%   This model handles the fluid equations of the splitting scheme and
%   setup a blackoil fluid model. The model is used in the fixed
%   stress splitting model.
%
% PARAMETERS:
%   G        - Grid
%   rock     - rock structure
%   fluid    - fluid structure
%   varargin -
%
% RETURNS:
%   class instance
%
% EXAMPLE:
%   run2DCase, runNorneExample
%
% SEE ALSO:
%   ThreePhaseBlackOilModel, MechFluidFixedStressSplitModel, MechFluidSplitModel

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
        function model = BlackOilFixedStressFluidModel(G, rock, fluid, varargin)
            model = model@ThreePhaseBlackOilModel(G, rock, fluid);
            model.disgas = true;
            model.vapoil = false;
            model = merge_options(model, varargin{:});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, ...
                                                        varargin)
            % Setup the equations for the fluid. The drivingForce contains
            % the volumetric changes computed from the mechanical equations.

            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  % Compatibility only

            opt = merge_options(opt, varargin{:});

            % Properties at current timestep
            [p, sW, sG, rs, rv, wellSol] = model.getProps(state, 'pressure', ...
                                                                 'water', ...
                                                                 'gas', ...
                                                                 'rs', 'rv', ...
                                                                 'wellSol');
            % Properties at previous timestep
            [p0, sW0, sG0, rs0, rv0] = model.getProps(state0, 'pressure', ...
                                                              'water', ...
                                                              'gas', 'rs', 'rv');

            %Initialization of primary variables ----------------------------------
            st  = model.getCellStatusVO(state,  1-sW-sG,   sW,  sG);
            st0 = model.getCellStatusVO(state0, 1-sW0-sG0, sW0, sG0);
            if model.disgas || model.vapoil
                % X is either Rs, Rv or Sg, depending on each cell's saturation status
                x = st{1}.*rs + st{2}.*rv + st{3}.*sG;
                gvar = 'x';
            else
                x = sG;
                gvar = 'sG';
            end

            [wellVars, wellVarNames, wellMap] = ...
                model.FacilityModel.getAllPrimaryVariables(wellSol);


            if ~opt.resOnly,
                % define primary varible x and initialize
                [p, sW, x, wellVars{:}] = ...
                    initVariablesADI(p, sW, x, wellVars{:});
            end

            fnew         = drivingForces.fixedStressTerms.new;
            mechTerm.new = fnew.sTerm + fnew.pTerm.*p;
            fold         = drivingForces.fixedStressTerms.old;
            mechTerm.old = fold.sTerm + fold.pTerm.*p0;

            otherDrivingForces = rmfield(drivingForces, 'fixedStressTerms');

            [eqs, names, types, state] = equationsBlackOilMech(state0, st0, ...
                                                              p, sW, x,  rs, ...
                                                              rv, st, wellVars, ...
                                                              state, model, ...
                                                              dt, mechTerm, ...
                                                              otherDrivingForces, ...
                                                              'iteration', ...
                                                              opt.iteration);

            primaryVars = {'pressure', 'sW', gvar, wellVarNames{:}};

            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end

        function forces = getValidDrivingForces(model)
            forces = getValidDrivingForces@ThreePhaseBlackOilModel(model);
            % divergence term
            % struct mechTerm.new and mechTerm.old
            forces.fixedStressTerms = [];
        end

        function fds = getAllVarsNames(model)
        % list of all the variable names that are used by this fluid model.
            fds = {'wellSol', 'pressure', 's', 'rs', 'rv'};
        end


    end

end
