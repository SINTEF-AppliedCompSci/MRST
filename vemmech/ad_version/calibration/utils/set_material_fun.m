function [funs, num_ctrls] = set_material_fun(G, props, matcond, u_offset)
%Undocumented Utility Function

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

   cells = (1:G.cells.num)';
   nlc = G.cells.nodePos(cells + 1) - G.cells.nodePos(cells);
   
   if isempty(matcond)
      funs.efun = @(arg) props.young;
      funs.nufun = @(arg) props.poisson;
      funs.loadfun = ...
          @(arg) SparseMultiArray(rldecode(props.rho(:), nlc), {'cn'}) * ...
                 SparseMultiArray(gravity(), {'d'});
      num_ctrls = [0, 0, 0];
      return
   end
   
   K_ctrls = find(matcond(:, 2) > 0); % controls for layer bulk modulus
   G_ctrls = find(matcond(:, 3) > 0); % controls for layer shear modulus
   R_ctrls = find(matcond(:, 4) > 0); % controls for layer density
   num_ctrls = [numel(K_ctrls), numel(G_ctrls), numel(R_ctrls)];
   
   layer_ixs = [matcond(:, 1); G.cartDims(3)+1];

   % determine the controls governing each cell in the model   
   cell_controls_K = distribute_layers(G, K_ctrls, layer_ixs);
   cell_controls_G = distribute_layers(G, G_ctrls, layer_ixs);
   cell_controls_R = distribute_layers(G, R_ctrls, layer_ixs);
   
   % compute the ranges within which we optimise
   K_sigma = matcond(K_ctrls, 2);
   G_sigma = matcond(G_ctrls, 3);
   R_sigma = matcond(R_ctrls, 4);
   
   K_ranges = [0 * K_sigma, K_sigma] + 1;
   G_ranges = [0 * G_sigma, G_sigma] + 1;
   R_ranges = [0 * R_sigma, R_sigma] + 1;   
   
   [K0, G0] = Enu2KG(props.young, props.poisson);
   
   k_rows = (1:G.cells.num)'; k_rows(cell_controls_K == 0) = [];
   k_mat = sparse(k_rows, cell_controls_K(cell_controls_K ~= 0), 1, G.cells.num, numel(K_ctrls));
   
   g_rows = (1:G.cells.num)'; g_rows(cell_controls_G == 0) = [];
   g_mat = sparse(g_rows, cell_controls_G(cell_controls_G ~= 0), 1, G.cells.num, numel(G_ctrls));
   
   rho_rows = (1:G.cells.num)'; rho_rows(cell_controls_R == 0) = [];
   rho_mat = sparse(rho_rows, cell_controls_R(cell_controls_R ~= 0), 1, G.cells.num, numel(R_ctrls));
   
   % entries are zero for cells with controls, 1 for cells without
   k_fixeds = ones(G.cells.num, 1) - k_mat * ones(num_ctrls(1), 1);
   g_fixeds = ones(G.cells.num, 1) - g_mat * ones(num_ctrls(2), 1);
   rho_fixeds = ones(G.cells.num, 1) - rho_mat * ones(num_ctrls(3), 1);
   
   % defining functions to specify E (Young's modulus) and nu (Poisson's parameter)
   efun = @(u) ...
     KG2E(K0 .* (zero_if_empty(k_mat * interp_range(u(1:num_ctrls(1)), ...
                                      K_ranges)) + k_fixeds), ...
          G0 .* (zero_if_empty(g_mat * interp_range(u(num_ctrls(1)+1:sum(num_ctrls(1:2))), ...
                                      G_ranges)) + g_fixeds));
   nufun = @(u)...
     KG2nu(K0 .* (zero_if_empty(k_mat * interp_range(u(1:num_ctrls(1)), ...
                                       K_ranges)) + k_fixeds), ...
          G0 .* (zero_if_empty(g_mat * interp_range(u(num_ctrls(1)+1:sum(num_ctrls(1:2))), ...
                                      G_ranges)) + g_fixeds));
   
   % defining load function (body forces)
   tmp = rldecode((1:G.cells.num)', nlc);
   
   loadfun = ...
       @(u) SparseMultiArray(...
           accumarray([(1:numel(tmp))', tmp], 1, [], [], [], true) * ...
           ( props.rho(:) .* ...
             (zero_if_empty(rho_mat * interp_range(u(end-num_ctrls(3)+1:end), R_ranges)) + ...
                     rho_fixeds) ), ...
           {'cn'}) * ...
       SparseMultiArray(gravity(), {'d'});
   
   % readjusting functions to ignore the controls that correspond to boundary
   % conditions
   funs.efun = @(u) efun(u(u_offset+1:end));
   funs.nufun = @(u) nufun(u(u_offset+1:end));
   funs.loadfun = @(u) loadfun(u(u_offset+1:end));
end

% ----------------------------------------------------------------------------
function res = zero_if_empty(var)
   if isempty(var)
      res = 0;
   else
      res = var;
   end
end

% ----------------------------------------------------------------------------
function res = interp_range(par, range)
   if isempty(value(par))
      res = [];
      return
   end
   par = par(:);
   res = (1-par) .* range(:,1) + par .* range(:,2);
end


% ----------------------------------------------------------------------------
function cell_layer = distribute_layers(G, controls, layer_ixs)
   
   layers = [layer_ixs(controls), layer_ixs(controls+1)-1];
   
   cell_layer = zeros(G.cells.num, 1);
   layer_size = prod(G.cartDims(1:2));
   
   for l = 1:size(layers, 1)

      start_ix = 1 + layer_size * (layers(l, 1) - 1);
      end_ix = layer_size * layers(l, 2);
      
      cell_layer(start_ix : end_ix) = l;
      
   end
   
   % @@will this always produce the intended result?
   cell_layer = cell_layer(G.cells.indexMap);
end

% ----------------------------------------------------------------------------
function [K, G] = Enu2KG(E, nu)
   K = E ./ (3 * (1 - 2 * nu));
   G = E ./ (2 + 2 * nu);
end
% ----------------------------------------------------------------------------

function E = KG2E(K, G)
   [E, ~] = KG2Enu(K, G);
end

% ----------------------------------------------------------------------------
function nu = KG2nu(K, G)
   [~, nu] = KG2Enu(K, G);
end

% ----------------------------------------------------------------------------
function [E, nu] = KG2Enu(K, G)

   E = 9 * K .* G ./ (3 * K + G);
   nu = (3 * K - 2 * G) ./ (2 * (3 * K + G));
end
