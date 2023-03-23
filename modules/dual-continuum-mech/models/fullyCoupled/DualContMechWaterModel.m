classdef DualContMechWaterModel < DualContMechFluidModel
%
% SYNOPSIS:
%   model = DualContMechWaterModel(G, rock, fluid, mech_problem, varargin)
%
% DESCRIPTION: 
%   Model for fully coupled mechanical single-phase (water) fluid simulation.
%
% PARAMETERS:
%   G     - Simulation grid.
%   rock  - Rock cell for fracture / matrix
%   fluid - Fluid cell for fracture / matrix
%
% RETURNS:
%   class instance
%
% EXAMPLE: consolidationTest
%
% SEE ALSO: DualContMechFluidModel
%
%{
Copyright 2009-2020 SINTEF ICT, Applied Mathematics.

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

    end

    methods
        function model = DualContMechWaterModel(G, rock, fluid, mech_problem, varargin)
            model = model@DualContMechFluidModel(G, rock, fluid, mech_problem, ...
                                                 varargin{:});              
        end

        function fluidModel = setupFluidModel(model)
            fluidModel = DualContWaterFluidModel(model.G, {model.rock, model.rock_matrix},...
                                                {model.fluid, model.fluid_matrix});
        end

        function [problem, state] = getEquations(model, state0, state, dt, ...
                                                        drivingForces, varargin)
            % Simulation options
            opt = struct('Verbose', mrstVerbose, ...
                         'reverseMode', false,...
                         'resOnly', false,...
                         'iteration', -1);  
            opt = merge_options(opt, varargin{:});

            % Shortcuts to sub-models
            fluidModel = model.fluidModel; 
            mechModel  = model.mechModel;  
                        
            % Extract variables and initialise as AD  properties at current timestep
            if ~opt.reverseMode
                [p, pm, wellSol, xd] = model.getProps(state, 'pressure', 'pressure_matrix',... 
                                                     'wellSol', 'xd');
                [wellVars, wellVarNames, ~] = fluidModel.FacilityModel.getAllPrimaryVariables(wellSol);                                 
                [p, pm, wellVars{:}, xd] = initVariablesADI(p, pm, wellVars{:}, xd);
            else
                error('Reverse mode AD currently not implemented for DC-mech module')
            end
            
            % well variables
            %[frac_index, mat_index] = fluidModel.findDCWells(wellSol);   
                        
            % Assemble linearised problem
            [w_eqs, w_eqsnames, w_eqstypes, state] = equationsDCWaterMech(state0, ...
                                                                          state, model, dt, ... 
                                                                          drivingForces, ...
                                                                         'iteration', ...
                                                                          opt.iteration);
                                                          
            [mech_eqs, mech_eqsnames, mech_eqstypes] = equationsDCPoroMechanics(xd, ...
                                                              mechModel, p, pm);                                                         
            eqs = horzcat(w_eqs, mech_eqs);
            names = {w_eqsnames{:}, mech_eqsnames{:}};
            types = {w_eqstypes{:}, mech_eqstypes{:}};
            primaryVars = {'pressure', 'pressure_matrix', wellVarNames{:}, 'xd'};
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        end       
        
        function mechTerm = computeStrainTerms(model, xd)
            % Method to compute strain terms (B_f:E and B_m:E) used in the
            % flow equations
            G = model.G;
            s = model.mechModel.operators;
            cM = model.mechModel.constitutive_coefficients_object;
            d = G.griddim;
            
            if any(diff(cM.B_f(:,1:d), 1, 2)) %anisotropic, would also work for B_m
                B_f = repmat(cM.B_f(:,1:d), 1, size(s.ovol_div,2)/d);
                B_f = B_f(:, ~s.isdirdofs);
                B_m = repmat(cM.B_m(:,1:d), 1, size(s.ovol_div,2)/d);
                B_m = B_m(:, ~s.isdirdofs);
            else %isotropic
                B_f = cM.B_f(:,1);
                B_m = cM.B_m(:,1); 
            end
            
            mechTerm.fracture = (s.div.*B_f)*xd./(G.cells.volumes);
            mechTerm.matrix = (s.div.*B_m)*xd./(G.cells.volumes);
            %mechTerm.fracture.old = (s.div.*B_f)*xd0./(G.cells.volumes);
            %mechTerm.matrix.old = (s.div.*B_m)*xd0./(G.cells.volumes);      
            % Note that the opmech.div returns the divergence integrated over cells. That is
            % why we divide by the cell's volumes
        end
       
    end
end
