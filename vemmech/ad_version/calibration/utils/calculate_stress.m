function stress = calculate_stress(u, dd, G, bcfun, efun, nufun, op)
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

   % make u an AD variable
   u = initVariablesADI(u);
      
   % computing (AD) values of dirichet degrees of freedoms (displacements at boundary)
   bc = bcfun(u);
   [u_bc, dirdofs] = bc.disp_bc.asVector({'d', 'n'}, [G.griddim, G.nodes.num]);
   skips = isnan(value(u_bc)); % flagged as 'free-moving', so not dirichlet
   dirdofs = setdiff(dirdofs, find(skips));
   
   % make dd into an AD variable, separate from u.  Restrict AD to those
   % elements of dd that are actual unknowns (i.e. not imposed displacements)
   dirdofs_ind = false(numel(dd), 1);
   dirdofs_ind(dirdofs) = true;
   dd_ad = initVariablesADI(dd(~dirdofs_ind));
   ndofs = numel(dd) - numel(dirdofs);
   dd_ad = sparse(find(~dirdofs_ind), 1:ndofs, 1, numel(dd), ndofs) * dd_ad;
   dd_ad(dirdofs) = dd(dirdofs); % these entries do not have any derivatives
   
   %% computing stress field as an AD variable only depending on dd
   tmp = op.WC' * (op.assemb' * dd_ad);

   stress_dd = op.D * tmp; 
   
   %% computing stress field as AD variable only depending on u
   
   % recompute operators as AD variables depending on u
   % @@ should be unneccessary since 'op' should already be right from when
   % 'dd' were calculated in the first place.  Hence, commented out.
   %[~, op] = VEM_assemble_AD(G, efun(u), nufun(u), 'extra', op);
   
   % make an AD version of 'dd' that only depends on 'u' (in dirichlet dofs)
   dd_u = repmat(u(1), numel(dd), 1) * 0; % establish correct AD structure
   dd_u(1:end) = dd;
   dd_u(dirdofs) = u_bc(dirdofs); % now, dirichlet entries are dependent on u
   
   tmp = op.WC' * (op.assemb' * dd_u); % tmp now dependent on u
   
   cellnum = size(op.WC, 2)/6; % @@ assumes 3D
   tmp = SparseMultiArray(tmp, ...
                      [repmat((1:6)', cellnum, 1), ...
                       rldecode((1:cellnum)', 6 * ones(cellnum, 1))], ...
                      {'j2', 'cell'});
   
   % stress field as AD variable only dependent on u
   stress_u = op.D_AD.product(tmp, true).contract('j', 'j2').asVector({'i', 'cell'});
   
   %% combining previous results into stress fields that depend both on u and dd

   stress = stress_u;
   stress.jac{2} = stress_dd.jac{1}; % now, 'stress' depends both on u and dd
   
   %% Putting stress field on the expected form to be returned

   % divide cross terms by two
   fac = ones(cellnum, 6);
   fac(:, 4:6) = 2;
   fac = reshape(fac', [], 1);
   stress = stress ./ fac;
   
   % swapping component 5 and 6 ( to be consistent with component order used
   % in 'calculateStressVEM')
   ix1 = 5:6:(cellnum*6);
   ix2 = 6:6:(cellnum*6);
   
   tmp = stress(ix1);
   stress(ix1) = stress(ix2);
   stress(ix2) = tmp;

end
