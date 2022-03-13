function [updata, report] = upAbsPerm(block, updata, varargin)
%Undocumented Utility Function

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
    'method',     'pressure', ...
    'dims',       1:3, ...
    'psolver',    'tpfa', ...
    'dp',         1*barsa ...
    );
opt = merge_options(opt, varargin{:});

if nargin==1
    updata = [];
end

if strcmpi(opt.method, 'pressure')
    f = @() upAbsPermPres(block, updata, 'dims', opt.dims, ...
        'psolver', opt.psolver, 'dp', opt.dp);
else
    f = @() upAbsPermAvg(block, updata, 'dims', opt.dims, ...
        'method', opt.method);
end

if nargout > 1
    [updata, report] = f();
else
    updata = f();
    report = [];
end

end
