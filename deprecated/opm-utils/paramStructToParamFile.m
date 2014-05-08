function paramStructToParamFile(s, filename)
   [fid, msg] = fopen(filename, 'wt');

   if fid < 0, error('Failed to open ''%s%'': %s', filename, msg); end

   % last line have to have lineshift
   prefix = '';
   writeParamStruct(fid, s, prefix);

   fclose(fid);
end

%--------------------------------------------------------------------------

function writeParamStruct(fid, s, prefix)
   fnames = fieldnames(s);

   for i = 1:numel(fnames),
      if isstruct(s.(fnames{i})),

         prefix_new = [prefix, '/', fnames{i}];
         writeParamStruct(fid, s.(fnames{i}), prefix_new);

      else

         val = s.(fnames{i});

         assert (sum([isnumeric(val), ...
                      islogical(val), ischar(val)]) == 1, ...
                 'Parameter value must be numeric, logical or string');

         if isnumeric(val),
            val = num2str(val);
         elseif islogical(val),
            t   = { 'false', 'true' };
            val = t{val + 1};
         end

         fprintf(fid, '%s/%s=%s\n', prefix, fnames{i}, val);
      end
   end
end
