function paper = paper_co2_optimalRatesPlacement()
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
paper = createPaperStruct('co2lab-optimalRatesPlacement', ...
                          'On obtaining optimal well rates and placement for CO2 storage', ...
                          'authors', 'R. Allen, H. M. Nilsen, O. Andersen, K.-A. Lie', ...
                          'published', 'Computational Geosciences, Vol. 21, Issue 5-6, pp. 1403-1422', ...
                          'year', 2017, ...
                          'modules', {'co2lab'}, ...
                          'doi', '10.1007/s10596-017-9631-6', ...
                          'url', 'https://link.springer.com/article/10.1007/s10596-017-9631-6');
end