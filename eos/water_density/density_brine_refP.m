function rho_b0 = density_brine_refP(rho_w0, t, mol_NaCl)
%Function that computes the density of brine at the reference pressure
% (p0 = 70Mpa). Equation 10 in Spivey et al., 2004

    T                   = t/100;

    % get the coefficients for pure water
    [ rho_w0, Ew, Fw ]  = coefficients_pure_water( t );

    % Table 3 in Spivey et al., 2004
    Dcm2_v              = [-7.925.*1e-5, -1.93.*1e-6, -3.42548.*1e-4, 0, 0];
    Dcm3_ov2_v          = [1.0998.*1e-3, -2.8755.*1e-3, -3.5819.*1e-3, -0.72877, 1.92016];
    Dcm1_v              = [-7.6402.*1e-3, 3.6963.*1e-2, 4.36083.*1e-2, -0.333661, 1.185685];
    Dcm1_ov2_v          = [3.746*.1e-4, -3.328.*1e-4, -3.346.*1e-4, 0, 0];


    % Compute the coefficients from Eq.(3) in Spivey et al., 2004
    Dcm2        = (Dcm2_v(1).*(T.^2) + Dcm2_v(2).*T + Dcm2_v(3)) ./ (Dcm2_v(4).*(T.^2) + Dcm2_v(5).*T +1);
    Dcm3_ov2    = (Dcm3_ov2_v(1).*(T.^2) + Dcm3_ov2_v(2).*T + Dcm3_ov2_v(3)) ./ (Dcm3_ov2_v(4).*(T.^2) + Dcm3_ov2_v(5).*T +1);
    Dcm1        = (Dcm1_v(1).*(T.^2) + Dcm1_v(2).*T + Dcm1_v(3)) ./ (Dcm1_v(4).*(T.^2) + Dcm1_v(5).*T +1);
    Dcm1_ov2    = (Dcm1_ov2_v(1).*(T.^2) + Dcm1_ov2_v(2).*T + Dcm1_ov2_v(3)) ./ (Dcm1_ov2_v(4).*(T.^2) + Dcm1_ov2_v(5).*T +1);


    % Compute density
    rho_b0      = rho_w0 + Dcm2.*(mol_NaCl.^2) + Dcm3_ov2.*(mol_NaCl.^(3/2)) + Dcm1.*mol_NaCl + Dcm1_ov2.*(mol_NaCl.^(1/2));

    tol = 1e-10;
    ix = abs(value(mol_NaCl)) < tol;
    rho_b0(ix) = rho_w0(ix);

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