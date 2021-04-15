function deck = readEDIT(fid, dirname, deck)
% Read edit

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

   [grd, miss_kw] = get_state(deck);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      switch kw
         case 'BOX'
            boxKeyword(fid);

         case 'ENDBOX'
            endboxKeyword;

         case 'MULTFLT'
            tmpl = { 'FaultName', '1.0', '1.0' };
            data = readDefaultedKW(fid, tmpl);  clear tmpl
            data(:, 2:end) = cellfun(@(s) sscanf(s, '%f'), ...
                                     data(:, 2:end),       ...
                                     'UniformOutput', false);

            if ~isfield(grd, kw), grd.(kw) = cell([0, 3]); end
            grd.(kw) = [grd.(kw); data];

         case {'DEPTH','PORV', 'TRANX', 'TRANY' , 'TRANZ' }
            grd = readGridBoxArray(grd, fid, kw, ...
                                   prod(deck.RUNSPEC.cartDims), NaN);

         case 'MULTPV'
            grd = readGridBoxArray(grd, fid, kw, ...
                                   prod(deck.RUNSPEC.cartDims), 1);

         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'}
            grd = applyOperator(grd, fid, kw);

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END'
            % Logical end of input deck.
            % Quite unusual (but nevertheless legal) in EDIT.
            %
            in_section = false;
            deck.GRID  = grd;

            % Restore default input box at end of section
            gridBox(defaultBox);

         case 'PROPS'
            % Read next section (always 'PROPS'.)
            in_section = false;

            deck = set_state(deck, grd, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readPROPS(fid, dirname, deck);

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, grd, miss_kw);

            deck = readEclipseIncludeFile(@readEDIT, fid, dirname, ...
                                          deck.RUNSPEC, deck);

            % Prepare for additional reading.
            [grd, miss_kw] = get_state(deck);

         otherwise
            if ischar(kw)
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, grd, miss_kw);
end

%--------------------------------------------------------------------------

function [grd, miss_kw] = get_state(deck)
   grd = deck.GRID;

   if ~isfield(deck.UnhandledKeywords, 'EDIT')
      miss_kw = {};
   else
      miss_kw = deck.UnhandledKeywords.EDIT;
   end
end

%--------------------------------------------------------------------------

function deck = set_state(deck, grd, miss_kw)
   deck.GRID                   = grd;
   deck.UnhandledKeywords.EDIT = unique(miss_kw);
end
