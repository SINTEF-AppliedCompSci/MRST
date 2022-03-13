function perm = upscalePerm(g, cg, rock, varargin)
% Compute upscaled permeabilites using flow based upscaling.
%
% SYNOPSIS:
%   Kc = upscalePerm(G, CG, rock)
%   Kc = upscalePerm(G, CG, rock, 'pn1', pv1, ...)
%
% PARAMETERS:
%   G, CG    - Underlying grid (G, see 'grid_structure') and coarse grid
%              model (CG) defined by function 'generateCoarseGrid'.
%
%   rock    - Rock data structure with valid field 'perm'.  The
%              permeability is assumed to be in measured in units of
%              metres squared (m^2).
%
%   'pn'/pv  - List of 'key'/value pairs defining optional parameters.
%
%                - method -- Discretization method used in computation.
%                         String: 'tpfa' (default) or 'mimetic'
%
%                - T   -- Transmissibility as defined by function
%                         'computeTrans' for the full fine-scale system.
%
%                - S   -- Mimetic linear system data structures as
%                         defined by function 'computeMimeticIP'.
%
%                - LinSolve --
%                         Handle to linear system solver software to which
%                         the fully assembled system of linear equations
%                         will be passed.  Assumed to support the syntax
%
%                             x = LinSolve(A, b)
%
%                         in order to solve a system Ax=b of linear eq.
%                         Default value: LinSolve = @mldivide (backslash).
%
%                - Verbose --
%                         Whether or not to emit progress reports during
%                         the assembly process.
%                         Logical.  Default value dependent upon global
%                         verbose settings of function 'mrstVerbose'.
%
% RETURNS:
%   Kc   - Upscaled permeability field
%
% REMARK:
%   The flow based upscaling applied here is based on the assumption that
%   the grid is Cartesian.
%
% SEE ALSO:
%   `extractSubgrid`, `computeTrans`, `computeMimeticIP`, `mrstVerbose`.

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

   try
      require mimetic incomp coarsegrid
   catch
      mrstModule add mimetic incomp coarsegrid
   end

   assert (isfield(g, 'cartDims'))

   opt = struct('Verbose',    mrstVerbose, ...
                'LinSolve',   @mldivide, ...
                'method',     'tpfa', ...
                'T', [], 'S', []);
   opt = merge_options(opt, varargin{:});

   if opt.Verbose,
      h = waitbar(0, 'Computing upscaled permeabilities...');
      fprintf('Computing upscaled permeabilities... '), tic
   end

   fluid    = initSingleFluid('mu' , 1 * Pascal*second, ...
                              'rho', 1 * kilogram/meter^3);
   dims     = find(g.cartDims);
   nBlocks  = cg.cells.num;
   perm     = zeros(cg.cells.num, numel(g.cartDims));

   switch lower(opt.method),
      case 'mimetic',
         if isempty(opt.S)
            s = computeMimeticIP(g, rock);
         else
            s = opt.S;
         end
         subS.ip   = s.ip; subS.type = s.type;
         ups = @(cells, hf) upscale_block_mimetic(...
            extractSubgrid(g, cells), s.BI(hf,hf), subS, dims, fluid, ...
            opt.LinSolve);
      case 'tpfa'
         if isempty(opt.T)
            T = computeTrans(g, rock);
         else
            T = opt.T;
         end
         ups = @(cells, hf) upscale_block_tpfa(...
            extractSubgrid(g, cells), T(hf), dims, fluid, opt.LinSolve);
   end

   map.sub_cells = @(b) find(cg.partition == b);
   map.hf        = @(c) mcolon(g.cells.facePos(  c  ), ...
                               g.cells.facePos(c + 1) - 1);

   % Compute upscaled permeability for each coarse block
   for b = 1 : nBlocks,
      cells = map.sub_cells(b);  % Fine-scale cells present in 'blk'
      hf    = map.hf(cells);     % Fine-scale half-faces in 'blk'

      perm(b,:) = ups(cells, hf);

      if opt.Verbose, waitbar(b / nBlocks, h), end
   end

   if opt.Verbose,
      toc, close(h)
   end
end


%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function k = upscale_block_mimetic(subG, BI, subS, dims, fluid, linsolve)
%
%  Compute perm in x-direction:
%
%                  v*n = 0
%              ______________
% l_boundary  |              | r_boundary
%             |              |
%     p1 = 1  |              | p2 = 0
%             |______________|
%                   v*n = 0
%
%             |--------------|
%                    L
%
%           grad p = (p2 - p1)/L
%           q1: flux over right boundary
%           kx* = (q1*L)/(area(r_boundary)*(p2 - p1))

   subS.BI = BI;
   k = zeros(1, numel(dims));
   rSol = initResSol(subG, 0.0);

   %{'LEFT', 'RIGHT', 'FRONT', 'BACK', 'BOTTOM', 'TOP'};
   tags = [ 1 2; 3 4; 5 6];

   bnd_f = any(subG.faces.neighbors==0,2);
   ind = bnd_f(subG.cells.faces(:,1));
   faceAndTag = subG.cells.faces(ind, :);

   for i = dims
      faces1 = faceAndTag(faceAndTag(:,2) == tags(i,1));
      faces2 = faceAndTag(faceAndTag(:,2) == tags(i,2));

      bc = addBC([], faces1, 'pressure', 1*barsa);
      bc = addBC(bc, faces2, 'pressure', 0);

      rSol = incompMimetic(rSol, subG, subS, fluid, 'bc', bc, ...
                           'LinSolve', linsolve);

      area = sum(subG.faces.areas(faces2,:));
      L    = abs(subG.faces.centroids(faces1(1), i) - ...
                 subG.faces.centroids(faces2(1), i));

      q    = abs(sum(rSol.flux(faces2)));
      k(i) = q * L / (1*barsa*area);
   end
end

function k = upscale_block_tpfa(subG, subT, dims, fluid, linsolve)
%
%  Compute perm in x-direction:
%
%                  v*n = 0
%              ______________
% l_boundary  |              | r_boundary
%             |              |
%     p1 = 1  |              | p2 = 0
%             |______________|
%                   v*n = 0
%
%             |--------------|
%                    L
%
%           grad p = (p2 - p1)/L
%           q1: flux over right boundary
%           kx* = (q1*L)/(area(r_boundary)*(p2 - p1))

   k = zeros(1, numel(dims));
   rSol = initResSol(subG, 0.0);

   %{'LEFT', 'RIGHT', 'FRONT', 'BACK', 'BOTTOM', 'TOP'};
   tags = [ 1 2; 3 4; 5 6];

   bnd_f = any(subG.faces.neighbors==0,2);
   ind = bnd_f(subG.cells.faces(:,1));
   faceAndTag = subG.cells.faces(ind, :);

   for i = dims
      faces1 = faceAndTag(faceAndTag(:,2) == tags(i,1));
      faces2 = faceAndTag(faceAndTag(:,2) == tags(i,2));

      bc = addBC([], faces1, 'pressure', 1*barsa);
      bc = addBC(bc, faces2, 'pressure', 0);

      rSol = incompTPFA(rSol, subG, subT, fluid, 'bc', bc, ...
                       'LinSolve', linsolve);

      area = sum(subG.faces.areas(faces2,:));
      L    = abs(subG.faces.centroids(faces1(1), i) - ...
                 subG.faces.centroids(faces2(1), i));

      q    = abs(sum(rSol.flux(faces2)));
      k(i) = q * L / (1*barsa*area);
   end
end
