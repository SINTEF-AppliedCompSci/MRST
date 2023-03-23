function [d_m, d_o, d_u, dpW] = getDispersionAnddpWMICP(model, state, poro)
%  Compute the disperison coefficients and "dpW" on the cell centers.

%{
Copyright 2021, NORCE Norwegian Research Centre AS, Computational 
Geosciences and Modeling. 

This file is part of the ad-micp module.

ad-micp is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ad-micp is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this file.  If not, see <http://www.gnu.org/licenses/>.
%}
  
  % Get model parameters
  alphaT = model.fluid.alphaT;
  alphaL = model.fluid.alphaL;
  Diff_m = model.fluid.diffm;
  Diff_o = model.fluid.diffo;
  Diff_u = model.fluid.diffu;
  muw = model.fluid.muw;
  % Compute current permeability
  K = model.fluid.K(poro);
  % "dpW" in the center of the cells
  v = faceFlux2cellVelocity(model.G, state.flux(:, 1));
  % Using the previous function to compute v results in half values of
  % velocity on the injection/free flow boundary. If the injetcion well 
  % is on the boundary, we correct the values
  if model.G.injectionwellonboundary == 1
      v(model.G.cellsinjectionwell, :) = ...
                                      2 * v(model.G.cellsinjectionwell, :);
  end
  v_norm = sqrt(sum(v .^ 2, 2));
  dpW = muw .* v_norm ./ K;
  % Dispersion terms for microbes, oxygen, and urea
  DD = (alphaL - alphaT) * ...
  [poro .* v(:, 1) .* v(:, 1) poro .* v(:, 1) .* v(:, 2) poro .* ...
       v(:, 1) .* v(:, 3) poro .* v(:, 2) .* v(:, 2) poro .* v(:, 2) .* ...
                     v(:, 3) poro .* v(:, 3) .* v(:, 3)] ./ (v_norm + eps);
  vm_tensor = (poro * Diff_m - alphaT * v_norm) * [1 0 0 1 0 1] + DD;
  vo_tensor = (poro * Diff_o - alphaT * v_norm) * [1 0 0 1 0 1] + DD;
  vu_tensor = (poro * Diff_u - alphaT * v_norm) * [1 0 0 1 0 1] + DD;
  % We use the harmonic average to compute the correspondic values on the
  % cell faces. We use the routine implemented in MRST to compute the
  % transmisibilities.
  Dispm = makeRock(model.G, vm_tensor, 1);
  Dispo = makeRock(model.G, vo_tensor, 1);
  Dispu = makeRock(model.G, vu_tensor, 1);
  d_m = getFaceTransmissibility(model.G, Dispm);
  d_o = getFaceTransmissibility(model.G, Dispo);
  d_u = getFaceTransmissibility(model.G, Dispu);
  % Half-trans -> trans and reduce to interior
  N = model.G.faces.neighbors;
  intInx = all(N ~= 0, 2);
  d_m = d_m(intInx);
  d_o = d_o(intInx);
  d_u = d_u(intInx);
end