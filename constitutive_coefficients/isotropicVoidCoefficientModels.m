classdef isotropicVoidCoefficientModels < isotropicStiffCoefficientModels
%
% SYNOPSIS:
%   model = isotropicVoidCoefficientModels(model)
%
% DESCRIPTION: 
%   Dervied class which instatiates constitutive coefficients 
%   using models from Khalili and Valliappan (1996) for which 
%   the underlying assumption is that the high-perm, low storage 
%   phase is all void space.
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
        function coefficients = isotropicVoidCoefficientModels(model)
            coefficients = coefficients@isotropicStiffCoefficientModels(model);
            
            % make sure the void space assumption is enforced
            phi_f = model.rock.poro./model.rock.vol_fraction;
            assert(all(phi_f == 1),...
                   'Intrinsic porosity for the void space coefficient models must equal 1')            
        end
        
        function [C_m, C_f] = caculateIntrinsicStiffnesses(coefficients, model)
        % Evaluate scalar (invariant) properties        
            C_m = Enu2C(model.mech.E_m, model.mech.nu_m, model.G);
            
            d = model.G.griddim;
            if d == 3 
                identity = bsxfun(@times, ones(model.G.cells.num, 1), [1, 1, 1, 0, 0, 0]);
            else 
                identity = bsxfun(@times, ones(model.G.cells.num, 1), [1, 1, 0]);
            end  
            
            K_m = (1/(d^2))*doubledot(identity, doubledot(identity, C_m, model.G), model.G);
            C_m = K_m;
            K_f = zeros(model.G.cells.num, 1);
            C_f = K_f;
        end      
    end
end

