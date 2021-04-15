function [sect, replaced] = replaceKeywords(sect, fid, keywords, nc)
% Replace keywords. Internal function

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

   replaced = false(size(keywords));

   kw = getEclipseKeyword(fid);
   in_section = ischar(kw);
   while in_section
      switch kw
         case {'ADD', 'COPY', 'EQUALS', 'MAXVALUE', ...
               'MINVALUE', 'MULTIPLY'}
            data = applyOperator(data, fid, kw);

         otherwise
            data = readGridBoxArray([], fid, kw, nc, 0.0);

            iskw = strcmpi(kw, keywords);

            if any(iskw)
               sect.(kw) = data.(kw);
               replaced(iskw) = true;
            end
      end

      kw = getEclipseKeyword(fid);
      in_section = ischar(kw);
   end
end
