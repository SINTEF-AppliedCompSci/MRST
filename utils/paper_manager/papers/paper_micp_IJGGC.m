function paper = paper_micp_IJGGC()
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
paper = createPaperStruct('micp_ijgcc', ...
    'Practical approaches to study microbially induced calcite precipitation at the field scale', ...
    'authors', 'D. Landa-Marbán, S. Tveit, K. Kumar and S. E. Gasda', ...
    'published', 'Int. J. Greenh. Gas Control 106, 103256', ...
    'year', 2021, ...
    'modules', {'ad-micp'}, ...
    'doi', '10.1016/j.ijggc.2021.103256', ...
    'url', 'https://doi.org/10.1016/j.ijggc.2021.103256');
end
