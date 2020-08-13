function stepinfo = readInfoStep(fname)
   [fid, msg] = fopen(fname, 'rt');
   if fid < 0
      error('Open:Failure', ...
            'Unable to Open INFOSTEP File ''%s'': %s', fname, msg);
   end
   clean = onCleanup(@() fclose(fid));
   vars = textscan(fgetl(fid), '%s');
   vars = regexprep(vars{1}, '\([^)]+\)', '');
   vals = textscan(fid, '%f');
   vals = num2cell(reshape(vals{1}, numel(vars), []) .', 1);
   stepinfo = cell2struct(vals, vars, 2);
end