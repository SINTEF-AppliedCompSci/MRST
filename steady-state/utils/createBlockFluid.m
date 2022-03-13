function BFluid = createBlockFluid(fluid, cells)
% Extracting the fluid for the current coarse cell only.

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

BFluid = fluid;

% Properties that take 'cellInx' name first
p = {'krW', 'krO', 'krOW', 'relPerm', 'pcOW', 'BW', 'bW', 'muW', ...
    'BO', 'bO', 'BOxmuO', 'muO', 'muWMult', 'ads', ...
    'pcOWInv', 'fracFlowInv'};
for i=1:numel(p)
    if isfield(fluid, p{i})
        BFluid.(p{i}) = @(x, varargin) fluid.(p{i})(x, ...
          'cellInx', subcells(cells, varargin{:}));
    end
end

end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function c = subcells(cells, varargin)
% This function is used to extract the correct cells from a fluid property.

opt = struct('cellInx', []);
opt = merge_options(opt, varargin{:});
if isempty(opt.cellInx)
   c = cells;
else
   c = cells(opt.cellInx);
end

end
