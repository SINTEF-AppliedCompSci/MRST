function paper = paper_MRST_book()
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
paper = createPaperStruct('mrst-book', ...
                          'An Introduction to Reservoir Simulation Using MATLAB; User Guide for the Matlab Reservoir Simulation Toolbox (MRST)', ...
                          'authors', 'K.-A. Lie', ...
                          'published', 'SINTEF ICT', ...
                          'year', 2015, ...
                          'modules', {'core', 'ad-core', 'ad-blackoil','book', 'coarsegrid','diagnostics','incomp','mimetic','mpfa','upscaling',}, ...
                          'url', 'http://www.sintef.no/projectweb/mrst/publications/', ...
                          'fileurl', 'http://www.sintef.no/contentassets/8af8db2e42614f7fb94fb0c68f5bc256/mrst-book-2015.pdf');
end