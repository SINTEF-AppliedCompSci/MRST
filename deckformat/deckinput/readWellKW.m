function w  = readWellKW(fid, w, kw)
%Read well definitions from an ECLIPSE Deck specification.
%
% SYNOPSIS:
%   w = readWellKW(fid, w, kw)
%
% PARAMETERS:
%   fid - Valid file identifier as obtained from FOPEN.
%
%   w   - Structure array containing current well configuration.  May be
%         empty if the well keyword is 'WELSPECS'.
%
%   kw  - Current well specification keyword.
%
% RETURNS:
%   w   - Updated structure array containing new data for keyword 'kw'.
%
% COMMENTS:
%   The currently recognized well keywords are:
%
%      'WELSPECS', 'COMPDAT', 'WCONHIST', 'WCONINJ', 'WCONINJE',
%      'WCONINJH', 'WCONPROD', 'GRUPTREE', 'WGRUPCON', 'GCONPROD',
%      'WPOLYMER', 'WELOPEN' and 'WELTARG'
%
%   Be advised that we do not support the complete feature set of these
%   keywords.
%
% SEE ALSO:
%   readEclipseDeck, processWells.

%{
Copyright 2009-2014 SINTEF ICT, Applied Mathematics.

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


   % Unless we're defining a new set of wells (or a well group hierarchy),
   % there had better be some previously defined wells.
   %
   assert (any(strcmp(kw, { 'GRUPTREE', 'WELSPECS' })) || ...
           ~isempty(w.WELSPECS), ...
          ['Well keyword ''%s'' encountered before any wells have ', ...
           'been declared using ''WELSPECS''.'], kw);

   switch kw,
      case 'WELSPECS', w = readWellSpec(fid, w);
      case 'COMPDAT' , w = readCompDat (fid, w);
      case 'WCONHIST', w = readWConHist(fid, w);
      case 'WCONINJ' , w = readWConInj (fid, w);
      case 'WCONINJE', w = readWConInje(fid, w);
      case 'WCONINJH', w = readWConInjh(fid, w);
      case 'WCONPROD', w = readWConProd(fid, w);
      case 'GRUPTREE', w = readGrupTree(fid, w);
      case 'GRUPNET' , w = readGrupNet (fid, w);
      case 'WGRUPCON', w = readWGrupCon(fid, w);
      case 'GCONINJE', w = readGConInje(fid, w);
      case 'GCONPROD', w = readGConProd(fid, w);
      case 'GECON'   , w = readGEcon   (fid, w);
      case 'WPOLYMER', w = readWPolymer(fid, w);
      case {'WELOPEN', 'WELLOPEN'},
         w = readWelOpen (fid, w);
      case {'WELTARG', 'WELLTARG'},
         w = readWelTarg(fid, w);
      otherwise
         fclose(fid);
         error(msgid('WellKW:Unsupported'), ...
               'Well keyword ''%s'' is not currently supported.', kw);
   end
end

%--------------------------------------------------------------------------
% Helpers follow.
%--------------------------------------------------------------------------

function w = readGrupTree(fid, w)
   %            1          2
   template = {'Default', 'Default'};
   gruptree = readDefaultedKW(fid, template);
   w.GRUPTREE = [w.GRUPTREE; gruptree];
end

%--------------------------------------------------------------------------

function w = readWellSpec(fid, w)
   %            1          2          3          4          5
   template = {'Default', 'Default', 'Default', 'Default', 'NaN', ...
      ...
      ...     % 6          7      8      9       10    11   12     13
               'Default', '0.0', 'STD', 'SHUT', 'NO', '0', 'SEG', '0'};
   numeric  = [3, 4, 5, 7, 11, 13];

   WellSpec = readDefaultedKW(fid, template);

   if ~isempty(WellSpec),
      WellSpec = toDouble(WellSpec, numeric);

      if ~isempty(w.WELSPECS),
         [tf, loc] = ismember(WellSpec(:,1), w.WELSPECS(:,1));
         ws = w.WELSPECS;
         if any(tf),
            ws(loc(tf),:) = WellSpec(tf,:);
         end
         WellSpec = [ws; WellSpec(~tf,:)];
      end
      w.WELSPECS = WellSpec;
   end
end

%--------------------------------------------------------------------------

function w = readCompDat(fid, w)
   %            1          2     3     4          5
   template = {'Default', '-1', '-1', 'Default', 'Default', ...
      ...
      ...     % 6       7     8       9      10
               'OPEN', '-1', '-1.0', '0.0', '-1.0', ...
      ...
      ...     % 11     12         13   14
               '0.0', 'Default', 'Z', '-1'};
   numeric  = [2:5, 7:11, 14];
   compdat  = toDouble(readDefaultedKW(fid, template), numeric);

   if isempty(w.COMPDAT),
      % No existing completion records.  Assign this data.
      w.COMPDAT = compdat;
   else
      % There are existing completions (typical case).  Need to decide
      % whether or not the 'compdat' records define completions of as-yet
      % unspecified wells (trivial), new completions of an existing
      % well (fairly easy) or modify existing completions of an existing
      % well (moderately involved).

      % 1) Identify wells for which existing records define completions.
      [ext.well, ie, ext.i] = unique(w.COMPDAT(:,1));                  %#ok

      % 2) Identify wells for which new records define completions.
      [new.well, in, new.i] = unique(  compdat(:,1));                  %#ok

      % Note: We do currently not support well lists ('*NAME') or well
      % templates ('NAME*') in this processing.
      if ~ all(cellfun(@isempty, regexp(new.well, '^\*|\*$', 'match'))),
         error('WLIST:Unsupp', ...
               'MRST does not support well lists or templates in COMPDAT');
      end

      [memb, loc] = ismember(new.well, ext.well);
      if ~ any(memb),
         % New records do not describe any existing completions.  Just
         % append the new records to the existing ones.
         w.COMPDAT = [ w.COMPDAT ; compdat ];
      else
         % New records affect existing wells.  Need to decide if these are
         % entirely new completions or if they modify some existing
         % completion data (e.g., to assign new well indices).  This is
         % ever so slightly involved, so defer processing to separate
         % helper.

         w = handleOverlapCompdat(w, compdat, memb, loc, ext, new);
      end
   end
end

%--------------------------------------------------------------------------

function w = readWConInj(fid, w)
   %            1          2          3       4          5      6
   template = {'Default', 'Default', 'OPEN', 'Default', 'inf', 'inf', ...
      ...
      ...     % 7      8       9     10     11     12
               '0.0', 'NONE', '-1', 'inf', '0.0', '0.0'};
   numeric  = [5:7, 9:numel(template)];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w = assignControlRecords(w, data, 'WCONINJ');
end

%--------------------------------------------------------------------------

function w = readWConInje(fid, w)
   %            1          2          3       4
   template = {'Default', 'Default', 'OPEN', 'Default', ...
      ...
      ...     % 5      6      7      8      9
               'inf', 'inf', 'NaN', 'inf', '0',         ...
      ...
      ...     % 10     11     12     13     14
               '0.0', '0.0', '0'  , '0'  , '0'};

   numeric  = 5 : numel(template);

   data = readDefaultedKW(fid, template);

   if ~isempty(data),
      data = toDouble(data, numeric);

      w = assignControlRecords(w, data, 'WCONINJE');
   end
end

%--------------------------------------------------------------------------

function w = readWConInjh(fid, w)
   %            1          2          3       4      5      6
   template = {'Default', 'Default', 'OPEN', '0.0', '0.0', '0.0', ...
      ...
      ...     % 7      8      9    10     11
               'inf', '0.0', '0', '0.0', '0.0'};
   numeric  = 4 : numel(template);

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w = assignControlRecords(w, data, 'WCONINJH');
end

%--------------------------------------------------------------------------

function w = readWConProd(fid, w)
   %            1          2       3          4      5      6
   template = {'Default', 'OPEN', 'Default', 'inf', 'inf', 'inf', ...
      ...
      ...     % 7      8      9      10     11   12
               'inf', 'inf', 'NaN', '0.0', '0', '0.0'};
   numeric  = 4 : numel(template);

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w = assignControlRecords(w, data, 'WCONPROD');
end

%--------------------------------------------------------------------------

function w = readWConHist(fid, w)
   %            1          2       3          4      5
   template = {'Default', 'OPEN', 'Default', '0.0', '0.0', ...
      ...
      ...     % 6      7    8          9      10     11
               '0.0', '0', 'Default', '0.0', '0.0', '0.0'};
   numeric  = [4 : 7, 9 : numel(template)];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w = assignControlRecords(w, data, 'WCONHIST');
end

%--------------------------------------------------------------------------

function w = readWelOpen(fid, w)
   %             1          2       3     4     5     6     7
   template = { 'Default', 'OPEN', '-1', '-1', '-1', '-1', '-1' };

   numeric = 3 : numel(template);

   data = readDefaultedKW(fid, template);

   assert (~ isempty(data), 'Internal error processing ''WELOPEN''.');

   data = toDouble(data, numeric);

   if ~ all([ data{:, 3:end} ] < 0),
      warning(['MRST does not support opening/closing individual ', ...
             'well connections through ''WELOPEN''.']);
   end

   for ikw = { 'WCONINJ', 'WCONINJE', 'WCONINJH' },
      if isfield(w, ikw{1}) && ...
            ~ (isempty(data) || isempty(w.(ikw{1}))),
         ix = ismember(data(:,1), w.(ikw{1})(:,1));

         injectors  = data(ix, :) .';
         data(ix,:) = [];

         for i = injectors,
            well = i{1};
            mode = i{2};

            ix = strcmp(well, w.(ikw{1})(:,1));
            w.(ikw{1}){ix, 3} = mode;
         end

         clear injectors ix
      end
   end

   for pkw = { 'WCONPROD', 'WCONHIST' },
      if isfield(w, pkw{1}) && ...
            ~ (isempty(data) || isempty(w.(pkw{1}))),
         ix = ismember(data(:,1), w.(pkw{1})(:,1));

         producers  = data(ix, :) .';
         data(ix,:) = [];

         for i = producers,
            well = i{1};
            mode = i{2};

            ix = strcmp(well, w.(pkw{1})(:,1));
            w.(pkw{1}){ix, 2} = mode;
         end

         clear producers ix
      end
   end

   assert (isempty(data), 'Internal error processing ''WELOPEN''.');
end

%--------------------------------------------------------------------------

function w = readGConInje(fid, w)
   %             1          2          3       4      5      6
   template = { 'Default', 'Default', 'NONE', 'Inf', 'Inf', 'Inf', ...
      ...
      ...     %  7      8      9      10   11     12     13
                'Inf', 'YES', 'NaN', '',  'NaN', 'NaN', 'Inf' };

   numeric  = [ 4 : 7 , 9 , 11 : 13 ];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   gconinje = w.GCONINJE;

   if ~ isempty(gconinje),
      [i, j] = blockDiagIndex(size(gconinje, 1), size(data, 1));

      m = strcmpi(gconinje(i, 1), data(j, 1)) & ...
          strcmpi(gconinje(i, 2), data(j, 2));

      if any(m),
         ij = [ i(m), j(m) ];

         assert (size(unique(ij, 'rows'), 1) == size(ij, 1), ...
                 'Unexpected match repetition in GECON');

         gconinje(ij(:,1), :) = data(ij(:,2), :);
         data    (ij(:,2), :) = [];
      end
   end

   w.GCONINJE = [ gconinje ; data ];
end

%--------------------------------------------------------------------------

function w = readGConProd(fid, w)
   %            1          2       3          4      5      6
   template = {'Default', 'NONE', 'Default', 'inf', 'inf', 'inf', ...
      ...
      ...     % 7       8      9      10  11      12      13
               'NONE', 'YES', 'inf', '', 'NONE', 'NONE', 'NONE', ...
      ...
      ...     % 14     15     16     17     18     19     20     21
               'inf', 'inf', 'inf', 'inf', 'inf', 'inf', 'inf', 'NONE'};
   numeric  = [3:6, 9, 14:20];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.GCONPROD = [w.GCONPROD; data];
end

%--------------------------------------------------------------------------

function w = readGEcon(fid, w)
   %             1          2      3      4      5      6
   template = { 'Default', '0.0', '0.0', '0.0', '0.0', '0.0', ...
      ...
      ...     %  7       8     9
                'NONE', 'NO', '0' };

   numeric  = [ 2 : 6, 9 ];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   gecon = w.GECON;

   if ~ isempty(gecon),
      [i, j] = blockDiagIndex(size(gecon, 1), size(data, 1));
      m      = strcmpi(gecon(i, 1), data(j, 1));

      if any(m),
         ij = [ i(m), j(m) ];

         assert (size(unique(ij, 'rows'), 1) == size(ij, 1), ...
                 'Unexpected match repetition in GECON');

         gecon(ij(:,1), :) = data(ij(:,2), :);
         data (ij(:,2), :) = [];
      end
   end

   w.GECON = [ gecon ; data ];
end

%--------------------------------------------------------------------------

function w = readWGrupCon(fid, w)
   %            1          2      3       4   5
   template = {'Default', 'YES', '-1.0', '', '1.0'};
   numeric  = [3,5];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WGRUPCON = [w.WGRUPCON; data];
end

%--------------------------------------------------------------------------

function w = readGrupNet(fid, w)
   %             1          2       3    4    5     6     7
   template = { 'Default', '-1.0', '0', '0', 'NO', 'NO', 'NONE' };
   numeric  = 2 : 4;

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   grpnet = w.GRUPNET;

   if ~ isempty(grpnet),
      [i, j] = blockDiagIndex(size(grpnet, 1), size(data, 1));
      m      = strcmpi(grpnet(i, 1), data(j, 1));

      if any(m),
         ij = [ i(m), j(m) ];

         grpnet(ij(:,1), :) = data(ij(:,2), :);
         data  (ij(:,2), :) = [];
      end
   end

   w.GRUPNET = [ grpnet ; data ];
end

%--------------------------------------------------------------------------

function w = readWPolymer(fid, w)
   %            1          2      3      4          5
   template = {'Default', '0.0', '0.0', 'Default', 'Default'};
   numeric  = [2, 3];

   data = readDefaultedKW(fid, template);
   data = toDouble(data, numeric);

   w.WPOLYMER = appendSpec(w.WPOLYMER, data, w.WELSPECS(:,1));
end

%--------------------------------------------------------------------------

function w = readWelTarg(fid, w)
   template = { 'Default', 'Default', 'NaN' };
   numeric  = numel(template);

   data = readDefaultedKW(fid, template);

   assert (~ isempty(data), 'Internal error processing ''WELTARG''.');

   data = toDouble(data, numeric);

   if ~ all(isfinite([ data{:,3} ])),
      error('MRST does not support defaulting Item 3 of ''WELTARG''');
   end

   for ikw = { 'WCONINJE', 'WCONINJH' },
      if isfield(w, ikw{1}) && ...
            ~ (isempty(data) || isempty(w.(ikw{1}))),
         inje = struct('RATE', 5, 'RESV', 6, 'BHP', 7);

         ix = ismember(data(:,1), w.(ikw{1})(:,1));
         injectors  = data(ix, :) .';
         data(ix,:) = [];

         fn = fieldnames(inje);
         for i = injectors,
            well = i{1};
            mode = i{2};

            if ~ any(strcmp(mode, fn)),
               warning('Inj:Targ:Unsupp', ...
                      ['Injection target mode ''%s'' unsupported', ...
                       '/ignored for well ''%s'' in ''WELTARG''.'], ...
                       mode, well);

               continue
            end

            ix = strcmp(well, w.(ikw{1})(:,1));
            w.(ikw{1}){ix, inje.(mode)} = i{3};
         end
      end
   end

   for pkw = { 'WCONPROD', 'WCONHIST' },
      if isfield(w, pkw{1}) && ...
            ~ (isempty(data) || isempty(w.(pkw{1}))),
         prod = struct('ORAT',  4, 'WRAT', 5, 'GRAT', 6, 'LRAT', 7, ...
                       'CRAT', 19, 'RESV', 8, 'BHP' , 9);

         ix = ismember(data(:,1), w.(pkw{1})(:,1));
         producers  = data(ix, :) .';
         data(ix,:) = [];

         fn = fieldnames(prod);
         for i = producers,
            well = i{1};
            mode = i{2};

            if ~ any(strcmp(mode, fn)),
               warning('Prod:Targ:Unsupp', ...
                      ['Production target mode ''%s'' unsupported', ...
                       '/ignored for well ''%s'' in ''WELTARG''.'], ...
                       mode, well);

               continue
            end

            ix = strcmp(well, w.(pkw{1})(:,1));
            w.(pkw{1}){ix, prod.(mode)} = i{3};
         end
      end
   end

   assert (isempty(data), 'Internal error processing ''WELTARG''.');
end

%--------------------------------------------------------------------------

function w = assignControlRecords(w, data, kw)
   assert (ischar(kw), 'Internal error');

   w.(kw) = appendSpec(w.(kw), data, w.WELSPECS(:,1));
   w      = remove(w, data(:,1), excludeSet(kw));
end

%--------------------------------------------------------------------------

function e = excludeSet(kw)
   e = { 'WCONINJ', 'WCONINJE', 'WCONINJH', 'WCONPROD', 'WCONHIST' };
   e = e(~ strcmpi(e, kw));
end

%--------------------------------------------------------------------------

function w = remove(w, name, exkw)
   for kw = reshape(exkw, 1, []),
      if isfield(w, kw{1}) && ~isempty(w.(kw{1})),
         w.(kw{1})(ismember(w.(kw{1})(:,1), name), :) = [];
      end
   end
end

%--------------------------------------------------------------------------

function T = toDouble(T, col)
   try
      T(:,col) = reshape(cellfun(@(v) sscanf(v, '%f'), T(:,col), ...
                                 'UniformOutput', false),        ...
                         [], numel(col));
   catch ME
      error('NonNumeric:Unexpected', ...
            'Unexpected non-numeric data encountered.');
   end
end

%--------------------------------------------------------------------------

function table = appendSpec(table, data, wells)
   matches = @(a,p) cellfun(@(c) ~isempty(c), regexp(a, p, 'match'));

   % Wildcard specs in new data?
   is_wc = matches(data(:,1), '\w+\*\s*$');
   for i = reshape(find(is_wc), 1, []),      % Row shape is essential here.
      patt = strrep(data{i,1}, '*', '.*');

      if ~isempty(table),
         % Remove any previous specs matching this 'patt'ern.
         table(matches(table(:,1), patt), :) = [];
      end

      match_wells = matches(wells, patt);
      nmatch      = sum(double(match_wells));
      if nmatch == 0,
         warning(msgid('Well:Unknown'), ...
                ['Well control wildcard ''%s'' does not match ', ...
                 'any previously defined wells.  Ignored.'], data{i,1});
      end
      append      = repmat(data(i,:), [nmatch, 1]);
      append(:,1) = wells(match_wells);

      table = [table; append];  %#ok       % Willfully ignore MLINT advice.
   end

   % Remove wildcard data.
   data(is_wc, :) = [];

   % Ensure that any remaining injection specs refer to known wells.
   if ~isempty(data) && ~all(ismember(data(:,1), wells)),
      unknown = data(~ismember(data(:,1), wells), 1);
      u = sprintf(' ''%s''', unknown{:});
      error(msgid('Well:Unknown'), ...
           ['Well control specified in undeclared wells:%s.\n', ...
            'Attempted coup d''état foiled.'], u);
   end

   % Exclude any previous records for these wells.  They were likely copied
   % during readSCHEDULE>defaultControl
   if ~ isempty(table),
      table(ismember(table(:,1), data(:,1)), :) = [];
   end

   % Complete keyword spec for this configuration.
   table = [table; data];
end

%--------------------------------------------------------------------------

function w = handleOverlapCompdat(w, compdat, memb, loc, ext, new)
   [ext.r, ext.pos] = compress(ext.i);  % Sparse rep. of exisisting records
   [new.r, new.pos] = compress(new.i);  % Sparse rep. of new records

   affected = false(size(ext.r));  % Existing records affected by overlap
   handled  = false(size(new.r));  % New records handled in overlap
   append   = cell([0, size(w.COMPDAT, 2)]);  % Accumulated result, overlap

   for o = reshape(find(memb), 1, []),
      wij = [ w.WELSPECS{strcmp(new.well{o}, w.WELSPECS(:,1)), 3:4} ];
      assert (numel(wij) == 2);

      % Overlapping wells
      or = ext.r(ext.pos(loc(o)) : ext.pos(loc(o) + 1) - 1);
      nr = new.r(new.pos(    o ) : new.pos(    o  + 1) - 1);

      assert (~ any(affected(or)), ...
              'This existing well (''%s'') already handled!?!', ...
              ext.well{loc(o)});
      assert (~ any(handled (nr)), ...
              'This new well (''%s'') already handled!?!', new.well{o});
      affected(or) = true;
      handled (nr) = true;

      ocd = expand(w.COMPDAT(or, :), wij);  % Existing COMPDAT for well
      ncd = expand(  compdat(nr, :), wij);  % New COMPDAT for well

      [i, j] = blockDiagIndex(size(ocd, 1), size(ncd, 1));
      oloc = reshape(vertcat(ocd{i, 2:5}), [], 4);
      nloc = reshape(vertcat(ncd{j, 2:5}), [], 4);

      k = all(oloc == nloc, 2);
      if any(k),
         % New records pertain to (some) existing records.  Update those.
         assert (size(unique([i(k), j(k)], 'rows'), 1) == sum(k), ...
                 'Multi-to-multi match not supported');

         ocd(i(k), :) = ncd(j(k), :);
      end

      append = [ append ; ocd ; ncd(unique(j(~k)), :) ];               %#ok
   end

   % New COMPDAT consists of those existing records that weren't affected,
   % followed by those new records that weren't handled in the overlap
   % region, followed by those overlapping records.  Note: This will likely
   % reorder the COMPDAT records, so the processing code in 'processWells'
   % *must* be robust in the face of scattered records.
   %
   w.COMPDAT = [ w.COMPDAT(~affected, :) ; compdat(~handled, :) ; append ];
end

%--------------------------------------------------------------------------

function [r, pos] = compress(range)
   assert (size(range, 2) == 1, 'range must be column vector');

   t   = sortrows([range, (1 : numel(range)) .']);
   r   = t(:,2);
   pos = cumsum([1 ; accumarray(t(:,1), 1) ]);
end

%--------------------------------------------------------------------------

function t = expand(t, ij)
   t(vertcat(t{:,2}) < 1, 2) = { ij(1) };   % I location
   t(vertcat(t{:,3}) < 1, 3) = { ij(2) };   % J location

   lo = vertcat(t{:,4});
   hi = vertcat(t{:,5});

   if any(hi > lo),
      % Single record accounts for more than one completion.  Expand to
      % 1-to-1 ratio of records-to-completions.

      t = rldecode(t, hi - lo + 1);
      i = num2cell(mcolon(lo, hi));

      t(:,4) = i;
      t(:,5) = i;
   end
end
