function [Kup, report] = upAbsPermPV(block, varargin)
% Simple pore-volume average of the absolute permeability
opt = struct(...
    'dims',       1:3 ...
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
pvSum = sum(pv);

% Compute average(s)
Kup   = nan(1,ndims);
if size(rock.perm, 2)==1
    % Isotropic permeability
    Kup(:) = sum(rock.perm(:,1).*pv)./pvSum; % same in all cells
else
    % Anisotropic permeability
    for i = 1:ndims % Loop over dimensions
        Kup(i) = sum(rock.perm(:,dims(i)).*pv)./pvSum;
    end
end

totalTime = toc(timeStart);
if wantReport
    report.dims     = dims;
    report.time     = totalTime;
end

end


