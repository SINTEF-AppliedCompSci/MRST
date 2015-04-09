function updata = upPcOW(block, updata, varargin)
% Upscale capillary pressure curves.
%
opt = struct(...
    'nvalues',  100 ...
    );
opt = merge_options(opt, varargin{:});

G     = block.G;
fluid = block.fluid;

if ~isfield(fluid, 'pcOW')
    updata.pcOW = [];
    return
end

% Find minimum and maximum pc in block
pc0   = fluid.pcOW( zeros(G.cells.num,1) );
pc1   = fluid.pcOW( ones(G.cells.num,1) );
pcmax = max(max(pc0), max(pc1));
pcmin = min(min(pc0), min(pc1));

% Pcow values to loop over
pcVec = linspace(pcmax, pcmin, opt.nvalues)';

pvSum  = sum(block.pv);
pcOWup = nan(opt.nvalues, 2);

% Loop over pc values and upscale
for ip = 1:opt.nvalues
    pcVal   = pcVec(ip);
    pc      = pcVal.*ones(G.cells.num, 1);
    pcOWInv = fluid.pcOWInv(pc); % saturation values in each cell
    sWup    = sum(pcOWInv.*block.pv) / pvSum;
    pcOWup(ip,:) = [sWup, pcVal];
end

% Add to structure
updata.pcOW = pcOWup;

end



