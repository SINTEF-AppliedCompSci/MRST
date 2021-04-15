function fluid = initSimpleScaledADIFluid(varargin)
% version of initSimpleScaledADIFluid with additional relperm scaling

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

if isa(varargin{1}, 'struct')
    % first input is fluid-structure -> update scaling
    fluid = varargin{1};
    if ~isfield(fluid, 'krW_base')
        opt = struct('swl', 0, 'swcr', 0, 'sowcr', 0, 'swu', 1);
        opt = merge_options(opt, varargin{2:end});
        fluid.krW_base  = fluid.krW;
        fluid.krO_base = fluid.krO;
    else
        opt = struct('swl', [], 'swcr', [], 'sowcr', [], 'swu', []);
        opt = merge_options(opt, varargin{2:end});
    end
    fn = fieldnames(opt);
    for k = 1:numel(fn)
        if ~isempty(opt.(fn{k}))
            fluid.(fn{k}) = opt.(fn{k});
        end
    end
else % initialize fluid
    opt = struct('swl', 0, 'swcr', 0, 'sowcr', 0, 'swu', 1);
    [opt, opfl] = merge_options(opt, varargin{:});
    fluid = initSimpleADIFluid(opfl{:});
    fluid.krW_base  = fluid.krW;
    fluid.krO_base = fluid.krO;
    fn = fieldnames(opt);
    for k = 1:numel(fn)
        fluid.(fn{k}) = opt.(fn{k});
    end
end

fluid.krW = @(s)fluid.krW_base(  scale_sw(s, fluid) );
fluid.krO = @(s)fluid.krO_base( scale_so(s, fluid) );
end

function s = scale_sw(s, f)
s = (s-f.swcr)./(f.swu-f.swcr);
s = max(0, min(1, s));
end

function s = scale_so(s, f)
s = (s-f.sowcr)./(1-f.swl-f.sowcr);
s = max(0, min(1, s));
end
