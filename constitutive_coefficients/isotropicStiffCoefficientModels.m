classdef isotropicStiffCoefficientModels < ConstitutiveCoefficients
%
% SYNOPSIS:
%   model = isotropicStiffCoefficientModels(model)
%
% DESCRIPTION: 
%   Dervied class which instatiates constitutive coefficients 
%   using models derived in Ashworth and Doster (TBR). Assumptions are
%   isotropy of the bulk material and the void space assumption. 
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
        function coefficients = isotropicStiffCoefficientModels(model)
            coefficients = coefficients@ConstitutiveCoefficients(model);                                 
        end
        
        function [C_m, C_f] = caculateIntrinsicStiffnesses(coefficients, model)
            % Evaluate scalar (invariant) properties 
            C_m = Enu2C(model.mech.E_m, model.mech.nu_m, model.G);
            C_f = Enu2C(model.mech.E_f, model.mech.nu_f, model.G);
            
            d = model.G.griddim;
            if d == 3 
                identity = bsxfun(@times, ones(model.G.cells.num, 1), [1, 1, 1, 0, 0, 0]);
            else 
                identity = bsxfun(@times, ones(model.G.cells.num, 1), [1, 1, 0]);
            end
            
            K_m = (1/(d^2))*doubledot(identity, doubledot(identity, C_m, model.G), model.G);
            C_m = K_m;
            K_f = (1/(d^2))*doubledot(identity, doubledot(identity, C_f, model.G), model.G);
            C_f = K_f;
        end 
        
        function [b_m, invn_m, b_f, invn_f] = calculateIntrinsicCoeffs(coefficients, model)
        % Calculate the intrinsic constitutive coefficients used in the 
        % effective parameter calculations. 
            % Intrinsic mechanical properties
            phi_f = model.rock.poro./model.rock.vol_fraction;
            phi_m = model.rock_matrix.poro./model.rock_matrix.vol_fraction;
            K_m = coefficients.iP.C_m;
            K_f = coefficients.iP.C_f; 
            
            % Intrinsic poroelastic properties
            b_m = ones(model.G.cells.num, 1) - K_m./(model.mech.K_s);
            invn_m = (b_m - phi_m)./(model.mech.K_s);
            b_f = ones(model.G.cells.num, 1) - K_f./(model.mech.K_s);
            invn_f = (b_f - phi_f)./(model.mech.K_s);
        end
        
        function [B_f, B_m, invN_f, invQ, invN_m] = caculateCoefficients(coefficients, model)
            % Calculate effective parameters
            [B_f, B_m, invN_f, invQ, invN_m] = ...
                calcIsotropicCoefficientModels(coefficients, model);
        end
    end
end

