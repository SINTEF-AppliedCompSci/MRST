function [updata, report] = upRelPermEPS(block, updata, ...
    method, varargin)
% Upscaling of relperm curves by end-point scaling (EPS)

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
    'dims',        1:3, ...  % Dimensions to upscale
    'dp',          1*barsa, ...  % Pressure drop
    'absmethod',   'pressure', ... % one-phase upscaling method
    'fullswof',    false, ...
    'nvalues',     [] ... % passed on to pv if fullswof=true
    );
opt = merge_options(opt, varargin{:});

wantReport  = nargout > 1;
timeStart   = tic;

% For now, we choose to limit the number of dimensions. This is just so we
% don't have to deal with it later. It's not really a limition here.
assert(numel(opt.dims)==1, 'Only supporting a single dim.');

% Assertions
assert(isa(block, 'GridBlock'), 'The block variable is not a GridBlock.');
assert(isfield(updata, 'perm'), 'One phase upscaling must be run first.');
assert(~isempty(method), 'Method must be given.');
assert(isprop(block, 'deck'), 'Need deck to do endpoint scaling.');

deck  = block.deck;
fluid = block.fluid;


%% Full SWOF upscaling if requested

% If requested, the SWOF tables are upscaled for every single block.
% Otherwise, there is only a single average SWOF table.
if opt.fullswof
    swopt  = {};
    if ~isempty(opt.nvalues)
        swopt = {'nvalues', opt.nvalues};
    end
    upswof = @(x) upRelPermPV(block, updata, 'porevolume', swopt{:});
    if wantReport
        [fullswofup, fullswofrep] = upswof();
    else
        fullswofup = upswof();
    end
end


%% SWIR and SOR upscaling

% Find values for each cell in block
swir  = cellfun(@(T) T(1,1), deck.PROPS.SWOF); % each region
swir  = swir(deck.REGIONS.SATNUM); % each cell in block
sor   = 1-cellfun(@(T) T(end,1), deck.PROPS.SWOF); % each region
sor   = sor(deck.REGIONS.SATNUM); % each cell in block
% krwm  = cellfun(@(T) T(end,2), block.deck.PROPS.SWOF); % each region
% krwm  = krwm(block.deck.REGIONS.SATNUM); % each cell in block

% Pore-volume average
pvTot = sum(block.pv);
pvf   = @(x) sum(x.*block.pv)./pvTot;
swirU = pvf(swir);
sorU  = pvf(sor);
% krwmU = pvf(krwm);


%% krW max upscaling

% Compute K*krW for maximum water saturation
krW  = fluid.krW(1-sor);
rockKkr = block.rock;
rockKkr.perm = bsxfun(@times, rockKkr.perm, krW);
block.rock = rockKkr;
Kup = updata.perm;

% Loop over the dimension
dims  = opt.dims;
ndims = length(dims);
krWmU = nan(1, ndims);
for id = 1:ndims
    d = dims(id); % Current dimension
    krKU = upAbsPerm(block, 'dims', d, 'dp', opt.dp, ...
                     'method', opt.absmethod);
    krWmU(id) = krKU / Kup(id); % Compute upscaled relperm value
end


%% pcOW(swir) upscaling

pcOW   = fluid.pcOW(swir);
pcOWmU = pvf(pcOW);


%% Wrap up

if opt.fullswof
    updata.pcOW = fullswofup.pcOW;
    updata.krO  = fullswofup.krO;
    updata.krW  = fullswofup.krW;
end

% Store upscaled data to structure
updata.swir    = swirU;
updata.sor     = sorU;
updata.krWmax  = krWmU;  % krW at 1-sor
updata.pcOWmax = pcOWmU; % pcOW at swir

% Report
if wantReport
    totalTime = toc(timeStart);
    report.method  = method;
    report.time    = totalTime;
    if opt.fullswof
        report.fullswof = fullswofrep;
    end
end


end




