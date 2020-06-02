classdef anisotropicCoefficientModels < ConstitutiveCoefficients
%
% SYNOPSIS:
%   model = anisotropicCoefficientModels(model)
%
% DESCRIPTION: 
%   Dervied class which instatiates constitutive coefficients 
%   using models derived in Ashworth and Doster (TBR). The only
%   assumption made is that the underlying phases and bulk material
%   are orthotropic (which can reduce to transfer isotropy) as necessary. 
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
% SEE ALSO: ConstitutiveCoefficients, DualContinuumMechanicModel
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
        function coefficients = anisotropicCoefficientModels(model)
            coefficients = coefficients@ConstitutiveCoefficients(model);                                 
        end
        
        function [C_m, C_f] = caculateIntrinsicStiffnesses(coefficients, model)
            C_m = Enu2C(model.mech.E_m, model.mech.nu_m, model.G);
            C_f = AnisoEnu2C(model.mech.E_f, model.mech.nu_f, model.mech.mu_f, model.G);
        end 
        
        function [b_m, invn_m, b_f, invn_f] = calculateIntrinsicCoeffs(coefficients, model)
        % Calculate the intrinsic constitutive coefficients used in the 
        % effective parameter calculations. 
            % Intrinsic mechanical properties
            phi_f = model.rock.poro./model.rock.vol_fraction;
            phi_m = model.rock_matrix.poro./model.rock_matrix.vol_fraction;
            C_m = coefficients.iP.C_m;
            C_f = coefficients.iP.C_f; 
            K_s = model.mech.K_s;          
            d = model.G.griddim;
            
            % The isotropic solid compliance S_si is defined as 
            % S_s:1 where 1 is the 2nd order identity tensor
            if d == 3
                S_si = bsxfun(@times, 1./(d.*K_s), [1, 1, 1, 0, 0, 0]);
                rank2 = @(x) bsxfun(@times, x, [1, 1, 1, 0, 0, 0]);
            else
                S_si = bsxfun(@times, 1./(d.*K_s), [1, 1, 0]);
                rank2 = @(x) bsxfun(@times, x, [1, 1, 0]);
            end
            
            % Intrinsic poroelastic properties
            b_m = rank2(ones(model.G.cells.num, 1)) - doubledot(C_m, S_si, model.G);
            invn_m = doubledot((b_m - phi_m), S_si, model.G);
            
            b_f = rank2(ones(model.G.cells.num, 1)) - doubledot(C_f, S_si, model.G);
            invn_f = doubledot((b_f - phi_f), S_si, model.G);    
        end
        
        function [B_f, B_m, invN_f, invQ, invN_m] = caculateCoefficients(coefficients, model)
            % Calculate effective parameters
            [B_f, B_m, invN_f, invQ, invN_m] = ...
                calcAnisotropicCoefficientModels(coefficients, model);
        end
    end
end

