function deck = readSUMMARY(fid, dirname, deck)
% Read summary

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


   [smry, miss_kw] = get_state(deck);

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section,
      switch kw,

         %-----------------------------------------------------------------
         % Sectioning keywords below.  Modifies flow of control.
         % Don't change unless absolutely required...
         %
         case 'END',
            % Logical end of input deck.
            %
            in_section   = false;
            deck.SUMMARY = smry;

         case 'SCHEDULE',
            % Read next section (always SCHEDULE)
            in_section = false;

            deck = set_state(deck, smry, miss_kw);

            % Restore default input box at end of section
            gridBox(defaultBox);

            deck = readSCHEDULE(fid, dirname, deck);

         case 'INCLUDE'
            % Handle 'INCLUDE' (recursion).
            deck = set_state(deck, smry, miss_kw);

            deck = readEclipseIncludeFile(@readSUMMARY, fid, dirname, ...
                                          deck.RUNSPEC, deck);

            % Prepare for additional reading.
            [smry, miss_kw] = get_state(deck);

         otherwise,
            if ischar(kw),
               miss_kw = [ miss_kw, { kw } ];  %#ok
            end
      end

      % Get next keyword.
      kw = getEclipseKeyword(fid);
      in_section = in_section && ischar(kw);
   end

   deck = set_state(deck, smry, miss_kw);
end

%--------------------------------------------------------------------------

function [smry, miss_kw] = get_state(deck)
   smry    = deck.SUMMARY;
   miss_kw = deck.UnhandledKeywords.SUMMARY;
end

%--------------------------------------------------------------------------

function deck = set_state(deck, smry, miss_kw)
   deck.SUMMARY                   = smry;
   deck.UnhandledKeywords.SUMMARY = unique(miss_kw);
end
