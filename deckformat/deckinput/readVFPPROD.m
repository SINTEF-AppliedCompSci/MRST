function [tid, tab] = readVFPPROD(fid)
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
