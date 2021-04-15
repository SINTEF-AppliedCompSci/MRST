function paper = paper_msmfe_slb()
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
paper = createPaperStruct('msmfe-slb', ...
    'Multiscale mimetic solvers for efficient streamline simulation of fractured reservoirs', ...
    'authors', 'J. R. Natvig, B. Skaflestad, F. Bratvedt, K. Bratvedt, K.-A. Lie, V. Laptev, and S. K. Khataniar', ...
    'published', 'SPE J., Vol. 16, No. 4, pp. 880-880', ...
    'year', 2011, ...
    'modules', {'msmfem'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/spe119132-PA.pdf', ...
    'doi', '10.2018/119132-PA', ...
    'url', 'https://www.onepetro.org/journal-paper/SPE-119132-PA');
end
