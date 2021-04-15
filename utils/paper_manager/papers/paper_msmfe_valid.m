function paper = paper_msmfe_valid()
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
paper = createPaperStruct('msmfe-valid', ...
    'Validation of the multiscale mixed finite-element method', ...
    'authors', 'M. Pal, S. Lamine, K.-A. Lie, and S. Krogstad', ...
    'published', 'Int. J. Numer. Meth. Fluids, Volume 77, Issue 4, pp. 206-223', ...
    'year', 2015, ...
    'modules', {'msmfe'}, ...
    'fileurl', 'http://folk.ntnu.no/andreas/papers/msmfem-validate.pdf', ...
    'doi', '10.1002/fld.3978', ...
    'url', 'http://onlinelibrary.wiley.com/doi/10.1002/fld.3978/abstract;jsessionid=52E2854CD309A48C538292E70A23570A.f01t03');
end
