function [data, pos] = insertInPackedData(pos, data, r, c)
%Insert c into row r of packed array (data, pos)
%
% SYNOPSIS:
%   [data, pos] = insertInPackedData(pos, data, r, c)
%
% PARAMETERS:
%   pos, data - Packed data array.
%   r         - Row indices of new data entries.
%   c         - New data entries.
%
% RETURNS:
%   data, pos - Updated packed data array.  Note reverse order compared to
%               input.
%
% EXAMPLE:
%   [G.cells.faces, G.cells.facePos] = ...
%       insertInPackedData(G.cells.facePos, G.cells.faces, ...
%                          cells, [faceIDs, hfTags])
%
%   [G.faces.nodes, G.faces.nodePos] = ...
%       insertInPackedData(G.faces.nodePos, G.faces.nodes, faces, newnodes)

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

   n   = numel(pos)-1;
   tmp = sortrows([r, c]);
   new = tmp(:,2:end);
   t   = accumarray(r, 1, [n, 1]);
   num = diff(pos);
   pos = cumsum([1;double(num)+t]);

   r   = unique(r);
   ix  = mcolon(pos(r)+num(r), pos(r+1)-1);
   newData = zeros(size(data, 1)+size(new, 1), size(data, 2));
   i    = (1:size(newData))';
   i(ix)=[];
   newData(i, :) = data;
   newData(ix,1:size(new,2)) = new;
   data = newData;
end
