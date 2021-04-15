function [tid, tab] = readVFPPROD(fid)
% Read VFP tables for producer

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

   r1 = { '0', '0.0', 'FLO', 'WFR', 'GFR', 'THP', 'ALQ', 'USYS', 'BHP' };

   data = readDefaultedRecord(fid, r1);

   data(1:2) = to_double(data(1:2));

   tid = data{1};
   tab = struct('tid'  , tid    , ...
                'depth', data{2}, ...
                'FLOID', data{3}, ...
                'WFRID', data{4}, ...
                'GFRID', data{5}, ...
                'THPID', data{6}, ...
                'ALQID', data{7}, ...
                'USYS' , data{8}, ...
                'QID'  , data{9});

   tab.FLO = readVector(fid, 'VFPPROD>FLO', inf);
   tab.THP = readVector(fid, 'VFPPROD>THP', inf);
   tab.WFR = readVector(fid, 'VFPPROD>WFR', inf);
   tab.GFR = readVector(fid, 'VFPPROD>GFR', inf);
   tab.ALQ = readVector(fid, 'VFPPROD>ALQ', inf);

   nflo = numel(tab.FLO);
   nthp = numel(tab.THP);
   nwfr = numel(tab.WFR);
   ngfr = numel(tab.GFR);
   nalq = numel(tab.ALQ);

   tab.Q = zeros([nflo, nthp, nwfr, ngfr, nalq]);
   nrec  = numel(tab.Q) / nflo;
   for rec = 1 : nrec,
      id = fscanf(fid, ' %d', 4);

      [thp, wfr, gfr, alq] = deal(id(1), id(2), id(3), id(4));

      tab.Q(:, thp, wfr, gfr, alq) = readVector(fid, 'VFPPROD>Q', nflo);
   end
end

%--------------------------------------------------------------------------

function v = to_double(v)
   convert = @(s) sscanf(regexprep(s, '[dD]', 'e'), '%f');

   if ischar(v),
      v = convert(v);
   else
      v = cellfun(convert, v, 'UniformOutput', false);
   end
end
