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
        upf = @(x) sum(x.*pv)/pvsum;
    case 'harmonic'
        upf = @(x) 1/( sum(pv./x)/pvsum );
    case 'geometric'
        % The geometric mean is defined as prod(x)^(1/numel(x)). However,
        % this is mathematically eqivalent to the expression below, and
        % the latter variant is preferred numerically.
        upf = @(x) exp( sum(pv.*log(x))/pvsum );
    otherwise
        error('Absolute permeability method unknown.');
end

% Compute average(s)
Kup   = nan(1,ndims);
if size(rock.perm, 2)==1
    % Isotropic permeability
    Kup(:) = upf( rock.perm(:,1) ); % same in all directions
else
    % Anisotropic permeability
    for i = 1:ndims % Loop over dimensions
        Kup(i) = upf( rock.perm(:,dims(i)) );
    end
end

totalTime = toc(timeStart);
if wantReport
    report.dims     = dims;
    report.time     = totalTime;
end

end


