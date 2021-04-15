function paper = paper_adjoint_msmfem()
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
paper = createPaperStruct('adjoint_msmfem', ...
                          'Adjoint Multiscale Mixed Finite Elements', ...
                          'authors', 'S. Krogstad, V.L Hauge, A. Gulbransen', ...
                          'published', 'SPE Journal, Vol. 16(1), pp. 162-171', ...
                          'year', 2011, ...
                          'doi', '10.2118/119112-PA', ...
                          'modules', {'adjoint','msmfem'}, ...
                          'url', 'http://dx.doi.org/10.2118/119112-PA');
end