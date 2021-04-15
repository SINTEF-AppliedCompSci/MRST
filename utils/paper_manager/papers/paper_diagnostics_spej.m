function paper = paper_diagnostics_spej()
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
paper = createPaperStruct('ms-diagnostics_spej', ...
                          'The application of flow diagnostics for reservoir management', ...
                          'authors', 'O. Moyner, S. Krogstad, and K.-A. Lie', ...
                          'published', 'SPE Journal, Vol. 20, No. 2, pp. 306-323', ...
                          'year', 2015, ...
                          'doi', '10.2118/171557-PA', ...
                          'modules', {'diagnostics'}, ...
                          'url', 'https://www.onepetro.org/journal-paper/SPE-171557-PA', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/diagnostics.pdf');
end