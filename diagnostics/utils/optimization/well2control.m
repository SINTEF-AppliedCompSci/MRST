function u = well2control(W, varargin)
% Produce control vector u from target-wells
% If scaling is supplied, u is scaled s.t. 0<=u<=1
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
opt  = merge_options(opt, varargin{:});
vals = vertcat(W(opt.targets).val);
if ~isempty(opt.scaling)
    lims = opt.scaling.boxLims;
    u = (vals - lims(:,1))./(lims(:,2)-lims(:,1));
    if any(or(u<-eps, u>1+eps))
        warning('Some well-values lie outside given box-constraints');
    end
else
    u = vals;
    warning('No scaling was given, setting control-values equal to target well-values');
end
end
