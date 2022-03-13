function bc = boundaryCellsSubGrid(G, c, varargin)
%Compute set of cells on the boundary of a subgrid
%
% SYNOPSIS:
%   bc = boundaryCellsSubGrid(G, c)
%
% PARAMETERS:
%   G - Grid structure.
%
%   c - List or logical mask of cells constituting a subgrid of 'G'.
%
%  'pn'/pv - List of 'key'/value pairs defining optional parameters.  The
%            supported options are:
%
%   neigbours - use other neibours than in grid
%
% RETURNS:
%   bc - List of cell indices that constitute the boundary cells of the
%        subgrid.

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

% $Date$
% $Revision$
opt = struct('neigbours',[]);
opt = merge_options(opt, varargin{:});

   if islogical(c),
      assert (numel(c) == G.cells.num, ...
              'Cell mask must be defined in all grid cells');
      c = find(c);
   end
  if(isempty(opt.neigbours))
    hf = mcolon(G.cells.facePos(c), G.cells.facePos(c + 1) - 1);
    li = zeros([G.cells.num + 1, 1]);   li(c + 1) = 1 : numel(c);
    N  = li(G.faces.neighbors(G.cells.faces(hf, 1), :) + 1);
    bc = c(accumarray(sum(N(sum(N > 0, 2) == 1, :), 2), 1) > 0);
  else
      N=opt.neigbours;
      p=false([G.cells.num+1,1]);
      p(c+1)=true;
      %Nred=false(size(N));
      Nred=p(N+1);
      orig_faces=sum(Nred,2)==1;
      N_tmp=double(N(orig_faces,:)).*double(Nred(orig_faces,:));
      bc=sum(N_tmp,2);
  end
   
end 