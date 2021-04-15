function tf = datasetHasCustomDownloadFcn(info)
%Predicate for whether or not a dataset provides a custom download function
%
% SYNOPSIS:
%   tf = datasetHasCustomDownloadFcn(info)
%
% PARAMETERS:
%   info - A dataset information structure as defined by function
%          datasetInfoStruct.
%
% RETURNS:
%   tf - Whether or not the dataset identified by 'info' does provide a
%        custom download function.

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

   tf = (numel(info.downloadFcn) == 1) && ...
         isa(info.downloadFcn{1}, 'function_handle') && ...
        (nargin(info.downloadFcn{1}) == 0);
end
