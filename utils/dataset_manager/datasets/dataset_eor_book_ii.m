function [info, present] = dataset_eor_book_ii()
% Info function for a set of water-EOR cases from the second MRST book

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
    [info, present] = datasetInfoStruct(...
        'name', 'EOR_Book_II', ...
        'description', 'Test cases for the ad-eor module', ...
        'fileurl', 'https://www.sintef.no/contentassets/124f261f170947a6bc51dd76aea66129/EOR_Book_II.zip', ...
        'hasGrid', true, ...
        'hasRock', true, ...
        'hasFluid', true, ...
        'filesize', 0.780 ...
         );
end
