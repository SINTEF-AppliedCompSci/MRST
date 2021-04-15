function paper = paper_agglom_flow()
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
paper = createPaperStruct('agglom-flow', ...
    'Flow-based coarsening for multiscale simulation of transport in porous media', ...
    'authors', 'V. L. Hauge, K.-A. Lie, and J. R. Natvig', ...
    'published', 'Comput. Geosci., Vol. 16, No. 2, pp. 391-408', ...
    'year', 2012, ...
    'modules', {'agglom'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/fb-grids.pdf', ...
    'doi', '10.1007/s10596-011-9230-x', ...
    'url', 'http://link.springer.com/article/10.1007%2Fs10596-011-9230-x');
end
