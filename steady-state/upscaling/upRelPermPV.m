function [updata, report] = upRelPermPV(block, updata, method, varargin)
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
    'nvalues',     20 ...
    );
opt = merge_options(opt, varargin{:});

wantReport  = nargout > 1;
timeStart   = tic;

ndims = 1; % all directions are equal, so we only use one dimension
nvals = opt.nvalues;
fluid = block.fluid;

assert(isa(block, 'GridBlock'));
assert(~isempty(method), 'Method must be set');

% Pore volume
pvTot = sum(block.pv);

% Allocate space
krW  = cell(1, ndims); krW(:) = {nan(nvals,2)};
krO  = cell(1, ndims); krO(:) = {nan(nvals,2)};
pcOW = nan(nvals, 2);

% Saturation values
satnum = block.deck.REGIONS.SATNUM;
sWmin  = cellfun(@(x) x(1,1), block.deck.PROPS.SWOF);
sWmax  = cellfun(@(x) x(end,1), block.deck.PROPS.SWOF);
sWmin  = sWmin(satnum);
sWmax  = sWmax(satnum);

values = linspace(0, 1, nvals)'; % sat.values unscaled

% Loop over the input values
for iv = 1:nvals
    val    = values(iv); % unscaled sW value
    sW     = val.*(sWmax-sWmin) + sWmin; % scaled
    sWup   = sum(sW.*block.pv)./pvTot;
    krOup  = sum(fluid.krOW(1-sW).*block.pv)./pvTot;
    krWup  = sum(fluid.krW(sW).*block.pv)./pvTot;
    pcOWup = sum(fluid.pcOW(sW).*block.pv)./pvTot;
    
    krO{1}(iv,:) = [sWup, krOup];
    krW{1}(iv,:) = [sWup, krWup];
    pcOW(iv,:)   = [sWup, pcOWup];
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
updata.krO  = krO;
updata.krW  = krW;
updata.pcOW = pcOW;

totalTime = toc(timeStart);
if wantReport
    report.method  = method;
    report.nvalues = nvals;
    report.time    = totalTime;
    report.valsOutsideRange = outside;
end


end




