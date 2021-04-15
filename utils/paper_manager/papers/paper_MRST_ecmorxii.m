function paper = paper_MRST_ecmorxii()
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
paper = createPaperStruct('mrst-ecmorxii', ...
                          'Discretisation on complex grids - Open source MATLAB implementation', ...
                          'authors', 'K.-A. Lie, S. Krogstad, I. S. Ligaarden, J. R. Natvig, H. M. Nilsen, and B. Skaflestad', ...
                          'published', 'Proceedings of ECMOR XII, Oxford, UK, 6-9 September 2010', ...
                          'year', 2010, ...
                          'modules', {'core', 'incomp', 'mimetic', 'mpfa', 'msmfem'}, ...
			  'doi', '10.3997/2214-4609.20145007', ...
                          'url', 'http://www.earthdoc.org/publication/publicationdetails/?publication=41326', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/ecmor-xii-mrst.pdf');
end