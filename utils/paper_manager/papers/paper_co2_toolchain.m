function paper = paper_co2_toolchain()
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
paper = createPaperStruct('co2lab-toolchain', ...
                          'An open-source toolchain for simulation and optimization of aquifer-wide CO2 storage', ...
                          'authors', 'O. Andersen, K.-A. Lie, H. M. Nilsen', ...
                          'published', 'Energy Procedia, Vol. 86, pp. 324-333', ...
                          'year', 2016, ...
                          'modules', {'co2lab'}, ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/Andersen-TCCS8.pdf', ...
                          'doi', '10.1016/j.egypro.2016.01.033', ...
                          'url', 'http://www.sciencedirect.com/science/article/pii/S1876610216000357');
end