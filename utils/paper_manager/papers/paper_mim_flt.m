function paper = paper_mim_flt()
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
paper = createPaperStruct('mim-flt', ...
    'Accurate modelling of faults by multipoint, mimetic, and mixed methods', ...
    'authors', 'H. M. Nilsen, K.-A. Lie, and J. R. Natvig', ...
    'published', ' SPE J., Vol. 17, No. 2, pp. 568-579', ...
    'year', 2012, ...
    'modules', {'mimetic'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/fault-multipliers.pdf', ...
    'doi', '10.2118/149690-PA', ...
    'url', 'https://www.onepetro.org/journal-paper/SPE-149690-PA');
end
