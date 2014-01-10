function pvtg = completePVTG(pvtg)
   for i = 1 : numel(pvtg),
      if any(diff(pvtg{i}.pos) == 1),
         error(['Scaled FVF copy not implemented for PVTG.\n', ...
                'Blame BSKA.']);
      end
   end
end
