classdef ConstitutiveCoefficients
%
% SYNOPSIS:
%   model = ConstitutiveCoefficients(G, model)
%
% DESCRIPTION: 
%   Base class for use in calculating constitutive coefficients belonging to the 
%   poroelastic dual-continuum constitutive model.  
%
% PARAMETERS:
%   G     - Grid structure
%   model - Takes the mechanical model
%
% RETURNS:
%   class instance
%
% EXAMPLE: 
%
% SEE ALSO: VoidSpaceCoefficients
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
        % Intrinsic parameters struct
        iP;
        
        % Effective constitutive coefficients
        B_f;
        B_m;
        invN_f;
        invQ;
        invN_m;
    end
    
    methods
        function coefficients = ConstitutiveCoefficients(model)
            % Check if we have the necessary fields required to calculate 
            % the intrinsic properties.                      
            assert(isfield(model.rock, 'poro'),...
                    'rock struct requires a fracture porosity field - "poro"')
            assert(isfield(model.rock_matrix, 'poro'),...
                    'rock_matrix struct requires a matrix porosity field - "poro"')
            assert(isfield(model.rock, 'vol_fraction'),...
                    'rock struct requires a fracture volume fraction field - "vol_fraction"')
            assert(isfield(model.rock_matrix, 'vol_fraction'),...
                    'rock_matrix struct requires a matrix volume fraction field - "vol_fraction"')
            assert(isfield(model.mech, 'K_s'),...
                    'mech struct requires a solid bulk modulus field - "K_s"')
                        
            % Calculate intrinsic stiffness - tensors in the anisotropic case,
            % moduli in the isotropic case. Assume C_m does not need to be calculated
            % a-priori since our underlying assumption is that of
            % homogeneous matrix material. 
            if ~isfield(model.mech, 'C_m') && ~isfield(model.mech, 'C_f')
                [coefficients.iP.C_m, coefficients.iP.C_f] = ...
                    caculateIntrinsicStiffnesses(coefficients, model);
            elseif ~isfield(model.mech, 'C_m') && isfield(model.mech, 'C_f')            
                [coefficients.iP.C_m, ~] = ...
                    caculateIntrinsicStiffnesses(coefficients, model);
            end
                
            % Calculate intrinsic coefficients
            [coefficients.iP.b_m, coefficients.iP.invn_m,...
                 coefficients.iP.b_f, coefficients.iP.invn_f]...
                = calculateIntrinsicCoeffs(coefficients, model); 
                                      
            % calculate effective coefficients
            [coefficients.B_f, coefficients.B_m, coefficients.invN_f, coefficients.invQ, coefficients.invN_m]...
                = caculateCoefficients(coefficients, model);
        end
        
        function [C_m, C_f] = caculateIntrinsicStiffnesses(coefficients, model)
            error('Base class function not meant for direct use.');
        end
        
        function [b_m, invn_m, b_f, invn_f] = caculateIntrinsicCoeffs(coefficients, model)
            error('Base class function not meant for direct use.');
        end
       
        function [B_f, B_m, invN_f, invQ, invN_m] = caculateCoefficients(coefficients, model)
            error('Base class function not meant for direct use.');
        end
    end
end

