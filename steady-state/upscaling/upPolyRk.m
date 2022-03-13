function [updata, report] = upPolyRk(block, updata, method, varargin)
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
    'nsat',        20, ... % Number of upscaled sat. values
    'npoly',       10, ... % Number of upscaled poly. values
    'verbose',  false  ...
    );
[opt, relpermOpt] = merge_options(opt, varargin{:});

% Also pass on nsat to relperm upscaling
relpermOpt = [relpermOpt {'nsat', opt.nsat, 'verbose', opt.verbose}];

G = block.G;
fluid = block.fluid;

% Values
ns = opt.nsat;
nc = opt.npoly;
svals = linspace(0.1, 1, ns)'; % unscaled values
cvals = linspace(0, fluid.cmax, nc)';

% Allocate
ndims = numel(updata.krW);
RkU = cell(1, ndims);
for d=1:3
    RkU{d}  = nan(ns, nc);
end
krWU = cell(1,ndims);

% One-phase upscaling data
up1.poro = updata.poro;
up1.perm = updata.perm;

dispif(opt.verbose, '\n\nStarting polymer Rk upscaling\n');
start = tic;

% Loop over polymer concentration values
for ic = 1:nc
    
    dispif(opt.verbose, '\nPolymer value %d of %d...\n', ic, nc);
    
    % Alter the fine-scale water relative permeability
    c  = cvals(ic).*ones(G.cells.num,1);
    Rk = 1 + (fluid.rrf-1).*( fluid.ads(c)./fluid.adsMax );
    block.fluid.krW = @(sW,varargin) fluid.krW(sW,varargin{:})./Rk;
    
    % Call the two-phase relative permeability upscaling
    [upd, report] = upRelPerm(block, up1, method, 'values', svals, ...
                              relpermOpt{:});
    
    % Loop over dimensions
    for d=1:ndims
        
        if ic == 1
            % c==0 is trivial to upscale, but we use the pure water krW
            % solution in the next iterations. This is to ensure that krW
            % and krW/Rk are upscaled for the same saturation values.
            krWU{d} = upd.krW{d};
            RkU{d}(:,ic) = 1;
        else
            % Compute upscaled Rk values
            RkU{d}(:,ic) = krWU{d}(:,2) ./ upd.krW{d}(:,2);
        end
        
    end
    
end

if opt.verbose
    sec = toc(start);
    min = floor(sec/60);
    sec = floor(mod(sec,60));
    timeStr = sprintf('%1.1f sec', sec);
    if min > 0, timeStr = [sprintf('%d min ', min) timeStr]; end
    dispif(opt.verbose, '\nCompleted upscaling of Rk in %s\n', timeStr);
end

sU = cell(1,ndims);
for i=1:ndims
    sU{i} = krWU{d}(:,1);
end

updata.Rk.val = RkU;
updata.Rk.s   = sU;
updata.Rk.c   = cvals;

end



