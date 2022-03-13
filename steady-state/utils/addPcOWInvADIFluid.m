function fluid = addPcOWInvADIFluid(fluid, deck, varargin)
% Add a fluid function for computing the inverse of the capillary pressure.
% This may be called later by using
%   sW = fluid.pcOWinv(pc)
% The input deck must have the property SWOF.

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

opt = struct(...
    'swir',  false  ... % include cravity in capillary part
	);
opt = merge_options(opt, varargin{:});

assert(isfield(deck, 'PROPS'), 'Invalid deck. Missing ''PROPS'' field.');
assert(isfield(deck.PROPS, 'SWOF'), 'Only property ''SWOF'' supported.');

reg   = handleRegions(deck, varargin{:});
swof  = deck.PROPS.SWOF;
fluid = assignPcOWInv(fluid, swof, reg);

end

% Remaining code is a slight modification of assignSWOF

function f = assignPcOWInv(f, swof, reg)
T = swof;
%T = cellfun(@(x)x(:,[4,1]), swof, 'UniformOutput', false);
haspcow = false;
for i=1:numel(T)
    if ~haspcow && any(T{i}(:,4) > 0)
        haspcow = true;
    end
    t = T{i}(:,[4,1]);
    if t(1,1)>t(2,1)
        t = flipud(t); % pcow values are decreasing, flip data
    end
    T{i} = extendTab(t);
end

% Only add pcOWInv if we have nonzero pcOW
if haspcow
    f.pcOWInv = @(pc, varargin)pcOWInv(pc, T, reg, varargin{:});
end

end


function v = pcOWInv(pc, T, reg, varargin)
    satinx = getRegMap(pc, reg.SATNUM, reg.SATINX, varargin{:});
    v = interpReg(T, pc, satinx);
end
