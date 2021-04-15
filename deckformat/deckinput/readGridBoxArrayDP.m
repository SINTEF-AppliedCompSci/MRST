function sect = readGridBoxArrayDP(sect, fid, kw, nc)
%Input grid array (global cell based) corresponding to current input box
%
% SYNOPSIS:
%   section = readGridBoxArray(section, fid, kw, nc)
%
% PARAMETERS:
%   section - Deck section data structure.  Expected to be one of GRID,
%             PROPS, REGIONS, or (in the future) SCHEDULE.
%
%   fid     - Valid file identifier as defined by function FOPEN.  The file
%             ponter FTELL(fid) is expected to be positioned directly
%             after the keyword that prompted reading of a grid array
%             (e.g., directly after the PERMX keyword).
%
%   kw      - Particular keyword (string) to which the next input data will
%             be assigned.  Input data will be extracted from the 'fid' for
%             all (global) cells corresponding to the current input box.
%
%   nc      - Number of cells in the global grid model.  Corresponds to
%             PROD(deck.RUNSPEC.DIMENS).  This number is used for
%             preallocating the output array ('kw') if it does not already
%             exist in the input deck section represented by 'section'.
%
% RETURNS:
%   section - Deck section data structure, modified such that new values
%             are assigned to the array 'fld' in the cells corresponding to
%             the current input box.
%
% SEE ALSO:
%   readGRID, readREGIONS, private/boxIndices, private/gridBox.

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


   fld = regexprep(kw, '\W', '_');

   if ~isfield(sect, fld), sect.(fld) = zeros([nc, 1]); end

   i = boxIndices;
    
   %rawRead= readVectorNew(fid, kw, nc);
    rawRead= readVectorNew(fid, kw, 'inf');
   

sect.(fld)(i,1) = rawRead(1:nc);
sect.(fld)(i,2) = rawRead((1+nc):end);
   
   
   
   
   
end