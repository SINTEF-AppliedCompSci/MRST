function partition = refineBlocks(partition, G, indicator, N_U, f_handle, varargin)
%Refine blocks for which indicator value exceeds given limit
%
% SYNOPSIS:
%   partition = refineBlocks(partition, G, indicator, N_U, ...
%                            @refineUniform, 'cartDims',[nx_c ny_c nz_c])
%   partition = refineBlocks(partition, G, indicator, N_U, ...
%                           @refineGreedy, 'neighbor_level', 1);
%
% DESCRIPTION:
%   This function refines too large blocks by use of the algorithm in the
%   function specified by f_handle.
%
% REQUIRED PARAMETERS:
%   partition - Partition vector
%
%   G         - Grid data structure discretising the reservoir model
%               (fine grid, geological model).
%
%   indicator - Cell-wise value of some measure/indicator function used for
%               deciding which blocks to refine.
%
%   indicator2- Cell-wise value of some measure/indicator function used for
%               deciding which neighboring block to merge into.
%
%   N_U       - Upper bound
%
% OPTIONAL PARAMETERS:
%
%   verbose  - Whether or not display number of blocks in the resulting
%              partition. Default value: verbose = false.
%
% RETURNS:
%   partition - Partition vector after refining.

% SEE ALSO:

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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


% Find, if there is one, the verbose option
i = find(strcmpi(varargin(1:2:end), 'verbose'));
if ~isempty(i),
   verbose = varargin{2*i};
   % If removing verbose from the varargin-list
   varargin = varargin([1:2*(i-1), 2*i+1:end]);
else
   verbose = false; % default
end

% Calling the actual refine-function
partition = f_handle(partition, G, indicator, N_U, varargin{:});

dispif(verbose, 'refineBlocks: %d blocks\n', max(partition));

end
