function [ax, ax_derivs] = compute_Ax_derivs(G, u, extra, E, nu, gamma, dirdofs)
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

   assert(G.griddim==3);
   N = G.cells.num;
   cpos = reshape(repmat(1:N, 6, 1), [], 1);
   rep6 = sparse((1:6*N)', cpos, 1, 6 * N, N); % repeat each vector entry 6 times

   nlc = diff(G.cells.nodePos);
   dim = 3;
   
   % repeat each entry '3n' times, where 'n' is the number of nodes for a given cell
   repX = sparse((1:dim*sum(nlc))', rldecode((1:N)', nlc * dim), 1); 
                                                                     
   D_dE = bsxfun(@rdivide, extra.D, rep6 * value(E));  %  d/dE (D)
   
   % we have D = Ds * Dm (where ds is scalar and dm is a matrix)
   Ds = E ./ ((1+nu) .* (1-2*nu)); 
   Dm = bsxfun(@rdivide, extra.D, rep6 * value(Ds));
   
   % computing d/dnu (D), which makes use of the definitions of Ds and Dm above
   Ds_dnu = E .* (1 + 4 * nu) ./ ( (1 + nu) .* (1 - 2 * nu) ).^2;
   
   Dm_dnu = spones(extra.D);
   Dm_diag = repmat([-1,-1,-1, -4,-4,-4], 1, G.cells.num);
   Dm_dnu(logical(eye(numel(Dm_diag)))) = Dm_diag;
   
   D_dnu = bsxfun(@times, Dm, rep6 * value(Ds_dnu)) + ...
           bsxfun(@times, Dm_dnu, rep6 * value(Ds));
   
   DNC = diag(extra.NC' * extra.NC);
   trDNC = sum(reshape(DNC, 6, []), 1)';
   c = gamma ./ trDNC .* G.cells.volumes; 
   trD = Ds .* 3 .* (3-5*value(nu)); 
   alpha_dE = c .* trD ./ value(E);
   alpha_dnu = - c .* value(E) .* 6 .* (value(nu)-1) .* (5 * value(nu) -1 ) ./ ...
                                    ( (value(nu) + 1) .* (1 - 2 * value(nu)) ).^2;
   S_dE = spdiags(repX * value(alpha_dE), 0, size(repX, 1), size(repX, 1)); 
   S_dnu = spdiags(repX * value(alpha_dnu), 0, size(repX, 1), size(repX, 1)); 
   
   % compute derivatives based on control variables
   ejac = []; nujac = [];
   if isa(E, 'ADI')
      ejac = [E.jac{:}];
   end
   if isa(nu, 'ADI')
      nujac = [nu.jac{:}];
   end
   if isempty(ejac)
      ejac = 0 * nujac;
   elseif isempty(nujac)
      nujac = 0 * ejac;
   end
   num_ders = size(ejac, 2);
      
   ImPP = extra.I - extra.PP;
   ax_derivs = sparse(numel(u), num_ders);
   for i = 1:num_ders
      
      % compute total derivative for each cell
      Kdu = ...
          extra.WC * ( ...
              bsxfun(@times, D_dE, rep6 * (ejac(:,i) .* G.cells.volumes)) + ...
              bsxfun(@times, D_dnu, rep6 * (nujac(:,i) .* G.cells.volumes))) * ...
          extra.WC' + ...          
          ImPP' * (...
              bsxfun(@times, S_dE, repX * ejac(:,i)) + ...
              bsxfun(@times, S_dnu, repX * nujac(:,i))) * ...
          ImPP;
      
      % assemble
      Kdu = extra.assemb * Kdu * extra.assemb';

      ax_derivs(:, i) = Kdu * u; % @@ how can we speed up here?
   end
   tmpd = dirdofs;
   dirdofs = false(size(ax_derivs, 1),1);
   dirdofs(tmpd) = true;
   ax_derivs = ax_derivs(~dirdofs, :);
   ax_derivs = full(ax_derivs);
   ax = extra.S * u;
   ax = ax(~dirdofs);
end
