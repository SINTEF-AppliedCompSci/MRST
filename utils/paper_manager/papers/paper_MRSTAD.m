function paper = paper_MRSTAD()
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
paper = createPaperStruct('mrst-ad', ...
                          'MRST-AD - an Open-Source Framework for Rapid Prototyping and Evaluation of Reservoir Simulation Problems', ...
                          'authors', 'S. Krogstad, K.-A. Lie, O. Moyner, H. M. Nilsen, X. Raynaud, and B. Skaflestad', ...
                          'published', 'Paper 173317-MS presented at the 2015 SPE Reservoir Simulation Symposium', ...
                          'year', 2015, ...
                          'modules', {'core', 'ad-core', 'ad-blackoil'}, ...
                          'url', 'https://www.onepetro.org/conference-paper/SPE-173317-MS', ...
                          'doi', '10.2118/173317-MS', ...
                          'fileurl', 'http://folk.ntnu.no/andreas/papers/rss15-mrst.pdf');
end