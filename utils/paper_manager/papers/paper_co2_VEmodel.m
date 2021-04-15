function paper = paper_co2_VEmodel()
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
paper = createPaperStruct('co2lab-VEmodel', ...
                          'Fully-implicit simulation of vertical-equilibrium models with hysteresis and capillary fringe', ...
                          'authors', 'H. M. Nilsen, K.-A. Lie, O. Andersen', ...
                          'published', 'Computational Geosciences, Vol. 20, No. 1, pp. 49-67', ...
                          'year', 2016, ...
                          'modules', {'co2lab'}, ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/co2lab-3.pdf', ...
                          'doi', '10.1007/s10596-015-9547-y', ...
                          'url', 'http://link.springer.com/article/10.1007%2Fs10596-015-9547-y');
end