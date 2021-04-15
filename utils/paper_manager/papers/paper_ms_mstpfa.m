function paper = paper_ms_mstpfa()
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
                          'published', 'Journal Computational Physics, Vol. 275, pp. 273-293', ...
                          'year', 2014, ...
                          'doi', '10.1016/j.jcp.2014.07.003', ...
                          'modules', {'msfvm', 'msrsb'}, ...
                          'url', 'http://www.sciencedirect.com/science/article/pii/S0021999114004835', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/mstpfa-paper.pdf');
end