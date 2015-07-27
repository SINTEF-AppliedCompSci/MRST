function CW = upWells(CG, rock, W, varargin)
% Upscale wells
opt = struct(...
    'LinSolve',   @mldivide ...
    );
opt = merge_options(opt, varargin{:});

% Compute transmissibility
s = setupSimComp(CG.parent, rock);
T = s.T_all;

% Perform upscaling
warning('off','upscaleTransNew:ZeroBoundary');
[~, ~, CW, ~] = upscaleTransNew(CG, T, ...
   'wells', {W}, 'bc_method', 'wells_simple', ...
   'LinSolve', opt.LinSolve );
warning('on','upscaleTransNew:ZeroBoundary');

end

