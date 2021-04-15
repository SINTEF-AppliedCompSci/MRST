function paper = paper_ms_msfvm()
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
paper = createPaperStruct('ms-msfvm', ...
                          'The multiscale finite-volume method on stratigraphic grids', ...
                          'authors', 'O. Moyner and K.-A. Lie', ...
                          'published', 'SPE Journal, Vol. 19, No. 5, pp. 816-831', ...
                          'year', 2014, ...
                          'doi', '10.2118/163649-PA', ...
                          'modules', {'msfvm'}, ...
                          'url', 'https://www.onepetro.org/journal-paper/SPE-163649-PA', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/SPE-163649-rev.pdf');
end