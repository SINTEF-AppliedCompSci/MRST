function smry = readEclipseSummary(files, smspec, open, readField)
%Read set of ECLIPSE summary files.
%
% SYNOPSIS:
%   summary = readEclipseSummary(files, smspec, open, readField)
%
% PARAMETERS:
%   files  - List, represented as a cell array of strings, of files from
%            which to extract summary data.
%
%   smspec - Summary specifiction as obtained from reading a '.FSMSPEC' or
%            '.SMSPEC' file.
%
%   open   - Read-open operation.  Must be compatible with fopen and
%            support the calling syntax
%
%                [fid, msg] = open(file)
%
%            the return values of which have the same semantics as the
%            equivalent 'fopen' return values.  In particular, 'fid' must
%            be a valid file identifier if the 'open' operation succeeds
%            and <0 otherwise.
%
%   readField -
%            Function handle offering the semantics of extracting a single
%            field (datum) from the summary data stream.  Assumed to
%            support the calling syntax
%
%                [name, data] = readField(fid, ...)
%
%            in which the 'fid' is a file stream obtained from an 'open'
%            call. The return value 'name' is assumed to be a single string
%            that names or identifies a particular datum, while 'data' is
%            the associated field data represented as a structure with the
%            fields '.type' and '.values'.
%
% NOTE:
%   If a summary data file contains multiple 'mini steps', then this data
%   will be concatenated in the appropriate fields.
%
%   If the field reader function does not match the summary data protocol,
%   e.g., if the data stream is from a formatted source while the field
%   reader expects unformatted (binary) data, then the behaviour of
%   function 'readEclipseSummary' is undefined.
%
% RETURNS:
%   summary - Summary data structure.  One field for each UNIQUE
%             smspec.KEYWORDS.values, and additional fields
%               - 'WNAME' -- Well or group names for active wells/groups.
%               - 'MINISTEP' --
%                            List of ministeps in set of summary files.
%
% SEE ALSO:
%   `private/readEclipseRestart`, `readEclipseOutputFileUnFmt`.

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


   smry = struct('RptTime', zeros([1, numel(files)]));
   k    = 1;
   for f = reshape(files, 1, []),
      [fid, msg] = open(f{1});
      if fid < 0, error([f{1}, ': ', msg]); end

      item = read_single_summary(fid, smspec, readField);

      fclose(fid);

      for i = reshape(fieldnames(item), 1, []),
         if any(strcmp(i{1}, {'UNITS', 'WNAME'})), continue, end

         field = i{1};
         if ~isfield(smry, field), smry.(field) = []; end

         smry.(field) = [smry.(field), item.(field)];
      end

      smry.RptTime(k) = numel(item.TIME);

      k = k + 1;
   end

   if ~isempty(smry),
      smry.RptTime = [1, cumsum(smry.RptTime)];
      if isfield(item, 'WNAME'), smry.WNAME = item.WNAME; end
      smry.UNITS = item.UNITS;
   end
end

%--------------------------------------------------------------------------

function smry = read_single_summary(fid, smspec, readField)
   ministep      = [];
   data          = struct();
   read_complete = false;
   % skip first 4
   fseek(fid, 4, 'cof');
   while ~read_complete;
      [name, field] = readField(fid);

      switch name,
         case 'SEQHDR',
            data.(name) = field;

         case 'MINISTEP',
            assert (numel(field.values) == 1, ...
                    'MINISTEP datum has %d values (expected 1).', ...
                    numel(field.values));
            stepno = field.values;

            [name, field] = readField(fid);
            assert (strcmp(name, 'PARAMS'), ...
                    'MINISTEP datum does not match expected field name.');

            if ~isfield(data, name), data.(name) = []; end

            data.(name) = [data.(name), field ];
            ministep    = [ministep   , stepno];  %#ok

         otherwise,
            if feof(fid),

               read_complete = true;

            elseif ~isempty(ferror(fid)),

               error('Input failure. System reports ''%s''.', ferror(fid));

            else

               error(msgid('Keyword:Unexpected'),                      ...
                    ['Unexpected keyword ''%s''.\nInput stream does ', ...
                     'not appear to be a valid summary file.'], name);

            end
      end
   end

   data.MINISTEP = ministep;
   smry = assign_keywords(data, smspec);
end

%--------------------------------------------------------------------------

function summary = assign_keywords(data, smspec)
   assert (isfield(smspec, 'KEYWORDS'), ...
           'Summary specification does not include ''KEYWORDS''.');

   [ukw, j, j] = unique(smspec.KEYWORDS.values);                       %#ok
   map         = sortrows([j, (1 : numel(j)) .']);
   pos         = cumsum([1; accumarray(map(:, 1), 1)]);
   kw_ix       = map(:, 2);

   summary = struct('MINISTEP', data.MINISTEP);

   if isfield(data, 'PARAMS'),
      vals = [ data.PARAMS.values ];
   else
      vals = zeros([numel(j), 1]);
   end

   for n = 1 : numel(ukw),
      field = regexprep(ukw{n}, '\W', '_');
      i     = kw_ix(pos(n) : pos(n + 1) - 1);

      % Unit of measurement (string) for this summary quantitiy.
      unit = unique(smspec.UNITS.values(i));  assert (numel(unit) == 1);

      summary.(field)       = vals(i, :);
      summary.UNITS.(field) = strtrim(unit{1});
   end

   % Extract well data only for defined (named) wells...
   matches = @(s,p) ~isempty(regexp(s, p, 'match', 'once'));
   wkw     = cellfun(@(s) matches(s, '^W'), ukw);

   if any(wkw),
      summary.WNAME = smspec.WGNAMES.values(j == find(wkw, 1, 'first'));
      is_act        = cellfun(@(s) matches(s, '[-\w]+'), summary.WNAME);
      ntotwell      = numel(is_act);

      for kw = reshape([ukw(wkw); 'WNAME'], 1, []),
         field = regexprep(kw{1}, {'+', '-'}, {'p', 'n'});
         field = regexprep(field, '\W', '_');

         if size(summary.(field), 1) == ntotwell,
            % Data items exists in result set for all declared wells.
            % Extract active only.  Otherwise, assume that the wanted data
            % items are explicitly stated in the run deck and do nothing...
            %
            summary.(field) = summary.(field)(is_act, :);
         end
      end
   end
end
