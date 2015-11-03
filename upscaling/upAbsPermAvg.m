function [Kup, report] = upAbsPermAvg(block, varargin)
% Power averaging of of the absolute permeability
opt = struct(...
    'dims',       1:3, ...
    'method',     'arithmetic' ...
    );
opt = merge_options(opt, varargin{:});

wantReport = nargout > 1;
timeStart = tic;

% Handle input
dims  = opt.dims; % Dimensions to upscale
ndims = length(dims);
G     = block.G;
rock  = block.rock;

% Pore-volume
pv    = rock.poro.*G.cells.volumes;
pvsum = sum(pv);

switch opt.method
    case 'arithmetic'
        upf = @(K,d) sum(K.*pv)/pvsum;
    case 'harmonic'
        upf = @(K,d) 1/( sum(pv./K)/pvsum );
    case 'geometric'
        % The geometric mean is defined as prod(x)^(1/numel(x)). However,
        % this is mathematically eqivalent to the expression below, and
        % the latter variant is preferred numerically.
        upf = @(K,d) exp( sum(pv.*log(K))/pvsum );
    case {'harmonic-arithmetic', 'arithmetic-harmonic'}
        upf = @(K,d) combinedAverage(pv, K, G.cartDims, d, opt.method);
    otherwise
        error('Absolute permeability method unknown.');
end

% Compute average(s)
Kup   = nan(1,ndims);

KupTMP = nan(1,ndims);

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

function K = combinedAverage(pv, perm, cartDims, d, method)
% Combined average of harmonic and arithmetic averages

% We do have an issue if some pore-volumes are zero. We set them to a
% tiny value instead.
pv(pv==0) = 10*eps;

PV = reshape(pv, cartDims);
K  = reshape(perm,  cartDims);

% The two other dimensions (not d)
dims    = 1:numel(cartDims);
odims   = [dims(1:d-1) dims(d+1:end)];

if strcmpi(method, 'harmonic-arithmetic')
    [K, PV] = harm(K, PV, d);
    [K, PV] = arit(K, PV, odims(1));
    [K, ~]  = arit(K, PV, odims(2));
elseif strcmpi(method, 'arithmetic-harmonic')
    [K, PV] = arit(K, PV, odims(1));
    [K, PV] = arit(K, PV, odims(2));
    [K, ~]  = harm(K, PV, d);
end
 
end

function [K, PV] = arit(K, PV, d)
% K and PV given as multidim matrix
pv = sum(PV, d);
K  = sum(PV.*K,d)./pv;  % arithmetic average
PV = pv;
end

function [K, PV] = harm(K, PV, d)
% K and PV given as multidim matrix
pv = sum(PV, d);
K  = 1./(sum(PV./K, d)./pv); % harmonic average
PV = pv;

% If the sum of the pore-volume is zero for a column in direction d,
% then the harmonic average of the permeability will be NAN. But as the
% pore-volume is zero, there will be no flow in this column, and we can
% set the permeability to zero.
K(isnan(K)) = 0;
end



