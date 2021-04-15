function paper = paper_co2_parameterUncertainty()
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
paper = createPaperStruct('co2lab-parameterUncertainty', ...
                          'Using simplified methods to explore the impact of parameter uncertainty on CO2 storage estimates with application to the Norwegian Continental Shelf', ...
                          'authors', 'R. Allen, H. M. Nilsen, K.-A. Lie, O. Møyner, O. Andersen', ...
                          'published', 'International Journal of Greenhouse Gas Control, Vol. 75, pp. 198-213', ...
                          'year', 2018, ...
                          'modules', {'co2lab'}, ...
                          'doi', '10.1016/j.ijggc.2018.05.017', ...
                          'url', 'https://www.sciencedirect.com/science/article/pii/S1750583617305753?via%3Dihub');
end