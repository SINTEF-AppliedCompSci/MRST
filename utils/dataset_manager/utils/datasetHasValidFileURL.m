function tf = datasetHasValidFileURL(info)
%Predicate for whether or not a dataset provides a valid file URL
%
% SYNOPSIS:
%   tf = datasetHasValidFileURL(info)
%
% PARAMETERS:
%   info - A dataset information structure as defined by function
%          datasetInfoStruct.
%
% RETURNS:
%   tf - Whether or not the dataset identified by 'info' does provide a
%        valid file URL from which to download the dataset.

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

   tf = (~ isempty(info.fileurl)) && ischar(info.fileurl);
end
