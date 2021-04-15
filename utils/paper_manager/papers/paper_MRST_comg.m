function paper = paper_MRST_comg()
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
paper = createPaperStruct('mrst-comg', ...
                          'Open source MATLAB implementation of consistent discretisations on complex grids', ...
                          'authors', 'K.-A. Lie, S. Krogstad, I. S. Ligaarden, J. R. Natvig, H. M. Nilsen, and B. Skaflestad', ...
                          'published', 'Computational Geosciences, , Vol. 16, No. 2, pp. 297-322', ...
                          'year', 2012, ...
                          'modules', {'core', 'mimetic', 'incomp', 'mpfa', 'msmfem'}, ...
			  'doi','10.1007/s10596-011-9244-4', ...
                          'url', 'http://link.springer.com/article/10.1007/s10596-011-9244-4', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/mrst-comg.pdf');
end