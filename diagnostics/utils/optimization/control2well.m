function W = control2well(u, W, varargin)
% Update val-fields of W from controls u
% If scaling is supplied, u is assumed to be scaled control s.t. 0<=u<=1
% according to scaling.boxLims

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

opt = struct('targets', (1:numel(W))', ...
             'scaling', []);
opt = merge_options(opt, varargin{:});
for k = 1:numel(opt.targets)
    wnr = opt.targets(k);
    if ~isempty(opt.scaling)
        bx = opt.scaling.boxLims(k,:);
        [umin, umax] = deal(bx(1), bx(2));
        W(wnr).val = u(k)*(umax-umin)+umin;
    else
        W(wnr).val = u(k);
        warning('No scaling was given, setting target well-values equal to control-values')
    end
end
end
