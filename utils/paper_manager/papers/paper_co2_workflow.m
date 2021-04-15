function paper = paper_co2_workflow()
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
paper = createPaperStruct('co2lab-workflow', ...
                          'A simulation workflow for large-scale CO2 storage in the Norwegian North Sea', ...
                          'authors', 'K.-A. Lie, H. M. Nilsen, O. Andersen, O. Moyner', ...
                          'published', 'Computational Geosciences, Vol. 20, No. 3, pp. 607-622', ...
                          'year', 2016, ...
                          'modules', {'co2lab'}, ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/co2lab-comg.pdf', ...
                          'doi', '10.1007/s10596-015-9487-6', ...
                          'url', 'http://link.springer.com/article/10.1007/s10596-015-9487-6'); 
end