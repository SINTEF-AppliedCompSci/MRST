function [B_f, B_m, invN_f, invQ, invN_m] = calcAnisotropicCoefficientModels(coefficients, model)
%
% SYNOPSIS:
%   function [B_f, B_m, invN_f, invQ, invN_m] = calcAnisotropicCoefficientModels(coefficients, model)
%
% DESCRIPTION: 
%   Computes anisotropic constitutive coefficients necessary for the
%   poroelastic dual-constitutive constitutive model using
%   the coefficient models presented in Ashworth and Doster 
%   (TBR)
%
% PARAMETERS:
%   coefficients - class instance
%   model 
%
% RETURNS:
%   B_f    - Effective Biot coefficient of the fracture
%   B_f    - Effective Biot coefficient of the matrix
%   invN_f - Effective solid Biot modulus of the fracture
%   invQ   - Inter-continuum pressure cross coupling coefficient
%   invN_m - Effective solid Biot modulus of the matrix
%
% EXAMPLE:
%
% SEE ALSO: calcIsotropicCoefficientModels, anisotropicCoefficientModels
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
C_dr = model.mech.C;
[v_m, v_f] = deal(model.rock_matrix.vol_fraction, model.rock.vol_fraction);
[C_m, b_m, invn_m, C_f, b_f, invn_f] = deal(coefficients.iP.C_m, coefficients.iP.b_m, ... 
                                            coefficients.iP.invn_m, coefficients.iP.C_f, ...
                                            coefficients.iP.b_f, coefficients.iP.invn_f);
invC_dif = computeInvC_dif(model.G, C_m, C_f);
                              
B_f = b_f - doubledot(b_f,...
            (doubledot((C_dr - C_f), invC_dif, model.G)), model.G); 
B_m = doubledot(b_m, (doubledot((C_dr - C_f), invC_dif, model.G)), model.G);
invN_f = doubledot(doubledot(b_f, invC_dif, model.G),...
                        (B_f - v_f.*b_f), model.G) + v_f.*invn_f;
invN_m = doubledot(doubledot(b_m, invC_dif, model.G),...
                        (v_m.*b_m - B_m), model.G) + v_m.*invn_m;
invQ = doubledot(doubledot(b_m, invC_dif, model.G),...
                        (v_f.*b_f - B_f), model.G);
                    
end

