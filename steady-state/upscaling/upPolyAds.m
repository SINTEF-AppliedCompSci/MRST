function updata = upPolyAds(block, updata, varargin)
% Upscale polymer adsorption isotherm using a simple average.
% 
% The adsorption isotherm which is a function of the local polymer solution
% concentration.
% 
opt = struct(...
    'npoly', 50 ...
    );
opt = merge_options(opt, varargin{:});

G = block.G;
f = block.fluid;
c = linspace(0, f.cmax, opt.npoly)';

% The adsorption is computed by a rock mass average
rm    = f.rhoR .* G.cells.volumes .* (1 - block.rock.poro);
rmtot = sum(rm);
ads   = nan(opt.npoly, 1);
for i=1:opt.npoly
   ads(i) = sum(rm .* f.ads( c(i).*ones(G.cells.num,1) ) ) / rmtot;
end

% Add to updata structure
updata.ads = [c ads];

end
