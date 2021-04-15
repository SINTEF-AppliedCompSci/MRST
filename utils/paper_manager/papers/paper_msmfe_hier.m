function paper = paper_msmfe_hier()
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
paper = createPaperStruct('msmfe-hier', ...
    'A hierarchical multiscale method for two-phase flow based upon mixed finite elements and nonuniform coarse grids', ...
    'authors', 'J. E. Aarnes, S. Krogstad and K.-A. Lie', ...
    'published', 'Multiscale Model. Simul., Vol. 5, No. 2, pp. 337-363', ...
    'year', 2006, ...
    'modules', {'msmfem'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/msmfem-grids.pdf', ...
    'doi', '10.1137/050634566', ...
    'url', 'http://epubs.siam.org/doi/abs/10.1137/050634566');
end
