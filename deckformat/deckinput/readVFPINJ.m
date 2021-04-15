function [tid, tab] = readVFPINJ(fid)
% Read VFP tables for injector

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

   r1 = { '0', '0.0', 'OIL', 'THP', 'USYS', 'BHP' };

   data = readDefaultedRecord(fid, r1);

   data(1:2) = to_double(data(1:2));

   tid = data{1};
   tab = struct('tid'  , tid    , ...
                'depth', data{2}, ...
                'FLOID', data{3}, ...
                'THPID', data{4}, ...
                'USYS' , data{5}, ...
                'BHPID', data{6});

   tab.FLO = readVector(fid, 'VFPINJ>FLO', inf);
   tab.THP = readVector(fid, 'VFPINJ>THP', inf);

   [nflo, nthp] = deal(numel(tab.FLO), numel(tab.THP));

   tab.BHP = zeros([nflo, nthp]);
   nrec    = nthp;
   for rec = 1 : nrec,
      thpid = fscanf(fid, ' %d', 1);

      tab.BHP(:, thpid) = readVector(fid, 'VFPINJ>BHP', nflo);
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
