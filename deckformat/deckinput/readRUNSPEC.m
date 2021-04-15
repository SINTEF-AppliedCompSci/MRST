function deck = readRUNSPEC(fid, dirname, deck)
% Read runspec

%Modiefied by VES for reading dual porosity eclipse decks


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

   [rspec, miss_kw] = get_state(deck);

   double_conv = @(c) sscanf(regexprep(c, '[dD]', 'e'), '%f');
   to_double = @(data) cellfun(double_conv, data);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      switch kw
         case 'AQUDIMS'
            tmpl = { '1', '1', '1', '36', '1', '1', '0', '0' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'COMPS'
            tmpl = {'0'};
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl

         case 'DIMENS'
            s = removeQuotes(readRecordString(fid));
            rspec.cartDims = reshape(sscanf(s, '%f', 3), 1, []);
            rspec.DIMENS   = rspec.cartDims;

            % Handle case of DIMENS *after* DUAL{PORO,PERM}
            if rspec.DUALPORO % if dual porosity then half the z cartdimes
                rspec.cartDims(3) = rspec.cartDims(3) / 2;
                rspec.DIMENS      = rspec.cartDims;
            end

            % Set default input box corresponding to entire model.
            defaultBox(rspec.DIMENS);

         case 'ENDSCALE'
            tmpl        = { 'NODIR', 'REVERS', '1', '20', '0' };
            data        = readDefaultedRecord(fid, tmpl);
            data(3:end) = num2cell(to_double(data(3:end)));   clear tmpl
            rspec.(kw)  = data;

         case 'EQLDIMS'
            tmpl = {'1', '100', '50', '1', '50'};
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl

         case 'EQLOPTS'
            data = readRecordString(fid);
            data = removeQuotes(tokenizeRecord(strtrim(data)));
            rspec.(kw)  = data;

         case 'GRIDOPTS'
            tmpl = { 'NO', '0', '0' };
            data = readDefaultedRecord(fid, tmpl);
            data(2:end) = cellfun(double_conv, data(2:end), ...
                                  'UniformOutput', false);
            rspec.(kw) = data;  clear tmpl data

         case 'MISCIBLE'
            tmpl = { '1', '20', 'NONE' };
            data = readDefaultedRecord(fid, tmpl);
            data(1:2) = cellfun(double_conv, data(1:2), ...
                                'UniformOutput', false);
            rspec.(kw) = data;  clear tmpl

         case 'PATHS'
            if ~isfield(rspec, kw), rspec.(kw) = cell([0, 2]); end

            tmpl       = { '|+|*ALIAS*|+|', '|+|*ACTUAL DIR*|+|' };
            rspec.(kw) = [ rspec.(kw) ; readDefaultedKW(fid, tmpl) ];

            clear tmpl

         case 'REGDIMS'
            tmpl = { '1', '1', '0', '0', '0', '1', '0', '0', '0' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear data tmpl

         case 'ROCKCOMP'
            tmpl = { 'REVERS', '1', 'NO', '', '0' };
            data = readDefaultedRecord(fid, tmpl);
            data([2,end]) = cellfun(double_conv, data([2,end]), ...
                                    'UniformOutput', false);
            rspec.(kw) = data;  clear tmpl data

         case 'SATOPTS'
            rspec.(kw) = readTextFlags(fid, rspec.(kw), kw);

         case 'START'
            s = readRecordString(fid);
            s = strrep(removeQuotes(s), 'JLY', 'JUL');
            rspec.START = datenum(s, 'dd mmm yyyy');

         case 'SMRYDIMS'
            rspec.(kw) = readVector(fid, kw, 1);  % Or just ignored...

         case 'TABDIMS'
            %        1    2     3     4     5     6     7    8
            tmpl = {'1', '1', '20', '20',  '1', '20', '20', '1', ...
               ...
               ... % 9    10     11    12   13   14   15    16
                    '1', 'NaN', '10', '1', '-1', '0', '0', 'NaN', ...
               ...
               ... % 17    18    19    20    21   22   23   24    25
                    '10', '10', '10', 'NaN', '5', '5', '5', '0', 'NaN', ...
                    };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl

         case 'TITLE'
            rspec.(kw) = strtrim(fgetl(fid));

         case 'UDADIMS'
            tmpl = { '0', '0', '100' };
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'UDQDIMS'
            tmpl = [{ '16', '16' }, repmat({ '0' }, [1, 8]), { 'N' }];
            data = readDefaultedRecord(fid, tmpl);
            data(1:end-1) = cellfun(double_conv, data(1:end-1), ...
                                       'UniformOutput', false);
            rspec.(kw) = data;   clear tmpl data

         case 'VFPIDIMS'
            % VFPIDIMS is a single record of up to three integer items.
            % Default values (=0) from E100.
            tmpl = repmat({ '0' }, [ 1, 3 ]);
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'VFPPDIMS'
            % VFPPDIMS is a single record of up to six integer items.
            % Default values (=0) from E100.
            tmpl = repmat({ '0' }, [ 1, 6 ]);
            data = readDefaultedRecord(fid, tmpl);
            rspec.(kw) = to_double(data);  clear tmpl data

         case 'WELLDIMS'
            %        1    2    3    4    5     6    7
            tmpl = {'0', '0', '0', '0', '5', '10', '5', ...
               ...
               ... % 8    9   10   11   12    13     14
                    '4', '3', '0', '1', '1', '10', '201'};
            data = readDefaultedRecord(fid, tmpl);  clear tmpl

            rspec.(kw) = to_double(data);
            if ~isfinite(rspec.(kw)(2))
               rspec.(kw)(2) = rspec.cartDims(3);
            end

         case 'WSEGDIMS'
            tmpl = { '0', '1', '1', '0' };
            data = readDefaultedRecord(fid, tmpl);

            rspec.(kw) = to_double(data);  clear data tmpl

         case {'NOGRAV', 'IMPES',        ...
               'METRIC', 'FIELD', 'LAB', 'PVT-M', ...
               'WATER' , 'OIL'  , 'GAS', 'SOLVENT', ...
               'DISGAS', 'VAPOIL', 'BLACKOIL', ...
               'POLYMER', 'SURFACT', 'BRINE', 'TEMP'}
            rspec.(regexprep(kw, '\W', '_')) = true;

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         % --------------- DP Keywords -------------- 
         case {'DUALPORO', 'DUALPERM', 'NODPPM'}
            rspec.(kw) = true;

            if strcmp(kw, 'DUALPERM')
                rspec.DUALPORO = true;
            end

            % Handle case of DUAL{PORO,PERM} *after* DIMENS
            if any(strcmp(kw, { 'DUALPORO', 'DUALPERM' })) && ...
                isfield(rspec, 'DIMENS')
                rspec.cartDims(3) = rspec.cartDims(3) / 2;
                rspec.DIMENS      = rspec.cartDims;
            end

            % Set default input box corresponding to entire model.
            defaultBox(rspec.DIMENS);

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END'
            % Logical end of input deck.
            % Quite unusual (but nevertheless legal) in RUNSPEC.
            %
            in_section   = false;
            deck.RUNSPEC = rspec;

            % Restore default input box at end of section
            gridBox(defaultBox);

         case 'GRID'
            % Read next section (always 'GRID').
            in_section = false;

            deck = set_state(deck, rspec, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readGRID(fid, dirname, deck);

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, rspec, miss_kw);

            deck = readEclipseIncludeFile(@readRUNSPEC, fid, dirname, ...
                                          deck.RUNSPEC, deck);

            % Prepare for additional reading.
            [rspec, miss_kw] = get_state(deck);

         otherwise
            if ischar(kw)
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, rspec, miss_kw);
end

%--------------------------------------------------------------------------

function [rspec, miss_kw] = get_state(deck)
   rspec   = initialise_dimension_arrays(deck.RUNSPEC);
   miss_kw = deck.UnhandledKeywords.RUNSPEC;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, rspec, miss_kw)
   deck.RUNSPEC                   = rspec;
   deck.UnhandledKeywords.RUNSPEC = unique(miss_kw);
end

%--------------------------------------------------------------------------

function rspec = initialise_dimension_arrays(rspec)
   if ~isfield(rspec, 'TABDIMS')
      rspec.TABDIMS = default_tabdims();
   end

   if ~isfield(rspec, 'WELLDIMS')
      rspec.WELLDIMS = default_welldims();
   end
end

%--------------------------------------------------------------------------

function tdims = default_tabdims()
   %        1  2   3   4  5   6   7  8
   tdims = [1, 1, 20, 20, 1, 20, 20, 1, ...
      ...
      ... % 9   10  11  12  13  14  15   16
            1, NaN, 10,  1, -1,  0,  0, NaN, ...
      ...
      ... % 17  18  19   20  21  22  23  24   25
            10, 10, 10, NaN,  5,  5,  5,  0, NaN];
end

%--------------------------------------------------------------------------

function wdims = default_welldims()
   %        1  2  3  4  5   6  7
   wdims = [0, 0, 0, 0, 5, 10, 5, ...
      ...
      ... % 8  9  10  11  12  13   14
            4, 3,  0,  1,  1, 10, 201];
end
