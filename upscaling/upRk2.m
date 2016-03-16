function [updata, report] = upRk2(block, updata, method, varargin)
opt = struct(...
    'nsat',        20,        ... % Number of upscaled sat. values
    'npoly',       10         ... % Number of upscaled poly. values
    );
[opt, relpermOpt] = merge_options(opt, varargin{:});

% Also pass on nsat to relperm upscaling
relpermOpt = [relpermOpt {'nsat', opt.nsat}];

G = block.G;
fluid = block.fluid;

ns = opt.nsat;
nc = opt.npoly;

svals = linspace(0, 1, ns)'; % Unscaled values
cvals = linspace(0, fluid.cmax, nc)';

swir = fluid.swir(fluid.satnum);
sor  = fluid.sor(fluid.satnum);

ndims = numel(updata.krW);
RkU = cell(1, ndims);
for d=1:3
    RkU{d}  = nan(ns, nc);
end

krWU = cell(1,ndims);

up1.poro = updata.poro;
up1.perm = updata.perm;

% Loop over polymer concentration values
for ic = 1:nc
    
    % Alter the fine-scale water relative permeability
    c  = cvals(ic).*ones(G.cells.num,1);
    Rk = 1 + (fluid.rrf-1).*( fluid.ads(c)./fluid.adsMax );
    block.fluid.krW = @(sW,varargin) fluid.krW(sW,varargin{:})./Rk;
    
    % Call the two-phase relative permeability upscaling
    [upd, report] = upRelPerm(block, up1, method, relpermOpt{:});
    
    % Loop over dimensions
    for d=1:ndims
        
        if ic == 1
            krWU{d} = upd.krW{d};
            RkU{d}(:,ic) = 1;
        else
            % Compute upscaled Rk values
            RkU{d}(:,ic) = krWU{d}(:,2) ./ upd.krW{d}(:,2);
            %RkU{d}(:,ic)
        end
        
    end
    
    
    disp('done');
    
end

sU = cell(1,ndims);
for i=1:ndims
    sU{i} = krWU{d}(:,1);
end

updata.Rk2.val = RkU;
updata.Rk2.s   = sU;
updata.Rk2.c   = cvals;

end



