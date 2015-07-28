function [updata, report] = upRelPermPV(block, updata, ...
    method, varargin)

opt = struct(...
    'nvalues',     20 ...
    );
opt = merge_options(opt, varargin{:});

wantReport  = nargout > 1;
timeStart   = tic;

ndims = 1; % all directions are equal, so we only use one dimension
nvals = opt.nvalues;
G     = block.G;
fluid = block.fluid;

assert(isa(block, 'GridBlock'));
assert(isfield(updata, 'perm'), 'One phase upscaling must be run first.');
assert(~isempty(method), 'Method must be set');

% Pore volume
pvTot = sum(block.pv);

% Allocate space
krW = cell(1, ndims);
krO = cell(1, ndims);

% Saturation values
sWmin = min(cellfun(@(x) x(1,1), block.deck.PROPS.SWOF));
sWmax = min(cellfun(@(x) x(end,1), block.deck.PROPS.SWOF));
values = linspace(sWmin, sWmax, nvals)';

% Loop over the input values
for iv = 1:nvals
    
    val   = values(iv);
    sW    = val.*ones(G.cells.num, 1);
    sWup  = sum(sW.*block.pv)./pvTot;
    krOup = sum(fluid.krOW(sW).*block.pv)./pvTot;
    krWup = sum(fluid.krW(sW).*block.pv)./pvTot;
    
    krO{1}(iv,:) = [sWup, krOup];
    krW{1}(iv,:) = [sWup, krWup];
end

% Check for upscaled values outside range. We simply force the values
% inside valid range.
outside = false;
for id = 1:ndims
    inx=krO{id}(:,2)>1; krO{id}(inx,2)=1; if any(inx),outside=true; end
    inx=krO{id}(:,2)<0; krO{id}(inx,2)=0; if any(inx),outside=true; end
    inx=krW{id}(:,2)>1; krW{id}(inx,2)=1; if any(inx),outside=true; end
    inx=krW{id}(:,2)<0; krW{id}(inx,2)=0; if any(inx),outside=true; end
end

% Store upscaled data to structure
updata.krO = krO;
updata.krW = krW;

totalTime = toc(timeStart);
if wantReport
    report.method  = method;
    report.nvalues = nvals;
    report.time    = totalTime;
    report.valsOutsideRange = outside;
end


end




