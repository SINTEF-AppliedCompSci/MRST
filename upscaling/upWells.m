function CW = upWells(CG, rock, W)
% Upscale wells


% Compute transmissibility
s = setupSimComp(CG.parent, rock);
T = s.T_all;

% Perform upscaling
warning('off','upscaleTransNew:ZeroBoundary');
[~, ~, CW, ~] = upscaleTransNew(CG, T, ...
  'wells', {W}, 'bc_method', 'wells_simple');
warning('on','upscaleTransNew:ZeroBoundary');

end


