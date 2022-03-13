function [updata, report] = upAbsPermAvg(block, updata, varargin)
% Power averaging of of the absolute permeability

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
    'dims',       1:3, ...
    'method',     'arithmetic' ...
    );
opt = merge_options(opt, varargin{:});

if nargin==1
    updata = [];
end

wantReport = nargout > 1;
timeStart = tic;

% Handle input
dims  = opt.dims; % Dimensions to upscale
ndims = length(dims);
G     = block.G;
rock  = block.rock;

% Volume
vol  = G.cells.volumes;
vsum = sum(vol);

switch opt.method
    case 'arithmetic'
        upf = @(K,d) sum(K.*vol)/vsum;
    case 'harmonic'
        upf = @(K,d) 1/( sum(vol./K)/vsum );
    case 'geometric'
        % The geometric mean is defined as prod(x)^(1/numel(x)). However,
        % this is mathematically eqivalent to the expression below, and
        % the latter variant is preferred numerically.
        upf = @(K,d) exp( sum(vol.*log(K))/vsum );
    case {'harmonic-arithmetic', 'arithmetic-harmonic'}
        upf = @(K,d) combinedAverage(vol, K, G.cartDims, d, opt.method);
    otherwise
        error('Absolute permeability method unknown.');
end

% Compute average(s)
Kup   = nan(1,ndims);

if size(rock.perm, 2)==1
    % Isotropic permeability
    if any(strcmpi(opt.method, {'harmonic-arithmetic', ...
            'arithmetic-harmonic'}))
        for i = 1:ndims % Loop over dimensions
            Kup(i) = upf( rock.perm, dims(i) );
        end
    else
        Kup(:) = upf( rock.perm(:,1) ); % same in all directions
    end
else
    % Anisotropic permeability
    for i = 1:ndims % Loop over dimensions
        Kup(i) = upf( rock.perm(:,dims(i)), dims(i) );
    end
end

updata.perm = Kup;

totalTime = toc(timeStart);
if wantReport
    report.method   = opt.method;
    report.dims     = dims;
    report.time     = totalTime;
end

end


%--------------------------------------------------------------------------
% HELPER FUNCTIONS
%--------------------------------------------------------------------------

function K = combinedAverage(vol, perm, cartDims, d, method)
% Combined average of harmonic and arithmetic averages

% We do have an issue if some volumes are zero. We set them to a tiny value
% instead.
vol(vol==0) = 10*eps;

V = reshape(vol, cartDims);
K = reshape(perm,  cartDims);

% The two other dimensions (not d)
dims    = 1:numel(cartDims);
odims   = [dims(1:d-1) dims(d+1:end)];

if strcmpi(method, 'harmonic-arithmetic')
    [K, V] = harm(K, V, d);
    [K, V] = arit(K, V, odims(1));
    [K, ~] = arit(K, V, odims(2));
elseif strcmpi(method, 'arithmetic-harmonic')
    [K, V] = arit(K, V, odims(1));
    [K, V] = arit(K, V, odims(2));
    [K, ~] = harm(K, V, d);
end
 
end

function [K, V] = arit(K, V, d)
% K and V given as multidim matrix
v = sum(V, d);
K = sum(V.*K,d)./v;  % arithmetic average
V = v;
end

function [K, V] = harm(K, V, d)
% K and PV given as multidim matrix
v = sum(V, d);
K = 1./(sum(V./K, d)./v); % harmonic average
V = v;

% If the sum of the pore-volume is zero for a column in direction d,
% then the harmonic average of the permeability will be NAN. But as the
% pore-volume is zero, there will be no flow in this column, and we can
% set the permeability to zero.
K(isnan(K)) = 0;
end
