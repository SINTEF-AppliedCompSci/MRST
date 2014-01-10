function [tid, tab] = readVFPINJ(fid)
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
