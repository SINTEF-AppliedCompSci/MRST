function deck = readREGIONS(fid, dirname, deck)
% Read regions

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

   [rgn, miss_kw] = get_state(deck);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      switch kw
         case 'BOX'
            boxKeyword(fid);

         case 'ENDBOX'
            endboxKeyword;

         case {'EQLNUM', 'FIPNUM' , 'IMBNUM' , ...
               'PVTNUM', 'SATNUM' , 'SURFNUM', ...
               'ENDNUM', 'ROCKNUM', 'FIPFAC' }

            nc = prod(deck.RUNSPEC.cartDims);
            if deck.RUNSPEC.DUALPORO
               rgn = readGridBoxArrayDP(rgn, fid, kw, nc);
            else
               rgn = readGridBoxArray(rgn, fid, kw, nc, 1);
            end

         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'}
            rgn = applyOperator(rgn, fid, kw);

         case {'ECHO', 'NOECHO'}
            kw = getEclipseKeyword(fid);
            continue;  % Ignore.  Not handled in MRST

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END'
            % Logical end of input deck.
            % Quite unusual (but nevertheless legal) in REGIONS.
            %
            in_section   = false;
            deck.REGIONS = rgn;

         case 'SOLUTION'
            % Read next section (always 'SOLUTION').
            in_section   = false;

            deck = set_state(deck, rgn, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readSOLUTION(fid, dirname, deck);

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, rgn, miss_kw);

            deck = readEclipseIncludeFile(@readREGIONS, fid, dirname, ...
                                          deck.RUNSPEC, deck);

            % Prepare for additional reading.
            [rgn, miss_kw] = get_state(deck);

         otherwise
            if ischar(kw)
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, rgn, miss_kw);
end

%--------------------------------------------------------------------------

function [rgn, miss_kw] = get_state(deck)
   rgn     = deck.REGIONS;
   miss_kw = deck.UnhandledKeywords.REGIONS;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, rgn, miss_kw)
   deck.REGIONS                   = rgn;
   deck.UnhandledKeywords.REGIONS = unique(miss_kw);
end
