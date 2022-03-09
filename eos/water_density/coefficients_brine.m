function [Eb, Fb] = coefficients_brine(Ew, Fw, mol_NaCl, t)
%Function that computes the coefficients for the brine equations (12)-(13)
%in Spivey et al., 2004

    T                = t/100;

    % Table 4 in Spivey et al., 2004
    Ecm_v            = [0, 0, 0.1353, 0, 0];
    Fcm3_ov_2_v      = [-1.409, -0.361, -0.2532, 0, 9.216];
    Fcm1_v           = [0, 5.614, 4.6782, -0.307, 2.6069];
    Fcm1_ov_2_v      = [-0.1127, 0.2047, -0.0452, 0, 0];

    % Compute coefficients Eq. (3) in Spivey et al., 2004
    Ecm             = (Ecm_v(1).*(T.^2) + Ecm_v(2).*T + Ecm_v(3)) ./ (Ecm_v(4).*(T.^2) + Ecm_v(5).*T +1);
    Fcm3_ov_2       = (Fcm3_ov_2_v(1).*(T.^2) + Fcm3_ov_2_v(2).*T + Fcm3_ov_2_v(3)) ./ (Fcm3_ov_2_v(4).*(T.^2) + Fcm3_ov_2_v(5).*T +1);
    Fcm1            = (Fcm1_v(1).*(T.^2) + Fcm1_v(2).*T + Fcm1_v(3)) ./ (Fcm1_v(4).*(T.^2) + Fcm1_v(5).*T +1);
    Fcm1_ov_2       = (Fcm1_ov_2_v(1).*(T.^2) + Fcm1_ov_2_v(2).*T + Fcm1_ov_2_v(3)) ./ (Fcm1_ov_2_v(4).*(T.^2) + Fcm1_ov_2_v(5).*T +1);

    % Compute Eb and Fb Eq. (12) and (13) in Spivey et al., 2004
    Eb              = Ew + Ecm.*mol_NaCl;
    Fb              = Fw + Fcm3_ov_2.*(mol_NaCl.^(3/2)) + Fcm1.*mol_NaCl + Fcm1_ov_2.*(mol_NaCl.^(1/2));

    tol = 1e-10;
    ix = abs(value(mol_NaCl)) < tol;
    Fb(ix) = Fw(ix);

end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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