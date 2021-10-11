function W = addAdjointWellFields(CG, W, varargin)
%
% SYNOPSIS:
%   W = addAdjointWellFields(CG, W, varargin)
%   W = addAdjointWellFields(CG, W, 'OverlapWell', 0, 'OverlapBlock', 4);
%
%
% DESCRIPTION:
%  Hack for the adjoint code:
%  Add additional fields to the Well-structure.
%  The added fields are coarseCells and optionally CS.overlap and
%  CS.wellOverlap.
%
%  Needed in updateWells and createSingleCellPseudoWells.

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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


opt = struct('OverlapWell', 0, 'OverlapBlock', 0);
opt = merge_options(opt, varargin{:});
overlapW = opt.OverlapWell;
overlapB = opt.OverlapBlock;

numWells = length(W);
for k = 1 : numWells,
   %{
   [trash, coarseCells] = find(CG.cells.subCells(W(k).cells,:));
   coarseCells       = unique(coarseCells);
   %}

   coarseCells = unique(CG.partition(W(k).cells));

   W(k).coarseCells  = coarseCells;

   W(k).CS.overlap     = overlapB;
   W(k).CS.wellOverlap = overlapW;
end
