function [d, vals, lims, flag] = extractSelectedDiagnostics(d, prop, tsel, wsel)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

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

%diagn = d.Data.diagnostics
%vals = zeros(size(D(1).tof, 1), ns);
%for k = 1:ns
[injIx, prodIx] = deal(wsel.injectorIx, wsel.producerIx);
if isempty(injIx)
   injIx = 1:numel(d.Data.diagnostics(tsel.ix(1)).D.inj);
end
if isempty(prodIx)
   prodIx = 1:numel(d.Data.diagnostics(tsel.ix(1)).D.prod);
end
flag = false; % true if 'vals' is a color a not a set of scalar values

switch prop
    case d.currentDiagnostics(1).name % TOF forward
        if isempty(d.currentDiagnostics(1).values)
            d = extractTOFForward(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(1).values;
        lims = [0 d.maxTOF/year];
        
    case d.currentDiagnostics(2).name % TOF backward
        if isempty(d.currentDiagnostics(2).values)
            d = extractTOFBackward(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(2).values;
        lims = [0 d.maxTOF/year];
        
    case d.currentDiagnostics(3).name % Residence time
        if isempty(d.currentDiagnostics(3).values)
            d = extractResidenceTime(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(3).values;
        lims = [0 2*d.maxTOF/year];
        
    case d.currentDiagnostics(4).name % Tracer forward
        if isempty(d.currentDiagnostics(4).values)
            d = extractTracerForward(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(4).values;
        lims = [0 1];
        
    case d.currentDiagnostics(5).name % Tracer backward
        if isempty(d.currentDiagnostics(5).values)
            d = extractTracerBackward(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(5).values;
        lims = [0 1];
        
    case d.currentDiagnostics(6).name % Tracer product
        if isempty(d.currentDiagnostics(6).values)
            d = extractTracerProduct(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(6).values;
        lims = [0 1];

    case d.currentDiagnostics(7).name % Sweep region
        if isempty(d.currentDiagnostics(7).values)
           d = extractSweepRegion(d, tsel.ix);
        end
        vals = d.currentDiagnostics(7).values;
        lims = [0 1];
        flag = true;

    case d.currentDiagnostics(8).name % Drainage region
        if isempty(d.currentDiagnostics(8).values)
           d = extractDrainageRegion(d, tsel.ix);
        end
        vals = d.currentDiagnostics(8).values;
        lims = [0 1];
        flag = true;
        
    case d.currentDiagnostics(9).name % First arrival forward
        if isempty(d.currentDiagnostics(9).values)
           d = extractFAForward(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(9).values;
        lims = [0 d.maxTOF/year];
        
    case d.currentDiagnostics(10).name % First arrival backward
        if isempty(d.currentDiagnostics(10).values)
           d = extractFABackward(d, tsel.ix, injIx, prodIx);
        end
        vals = d.currentDiagnostics(10).values;
        lims = [0 d.maxTOF/year];
        
    otherwise
        warning(['Diagnostics property ', prop, ' is not known or badly spelled....'])
end
end

function d = extractTracerForward(d, ts, injIx, ~)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
v = zeros(d.G.cells.num, nsteps);
for k = 1:nsteps
    v(:, k) = sum(diagn(ts(k)).D.itracer(:, injIx), 2);
end
d.currentDiagnostics(4).values = v;
end

function d = extractTracerBackward(d, ts,  ~, prodIx)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
v = zeros(d.G.cells.num, nsteps);
for k = 1:nsteps
    v(:, k) = sum(diagn(ts(k)).D.ptracer(:, prodIx), 2);
end
d.currentDiagnostics(5).values = v;
end

function d = extractTOFForward(d, ts, injIx, prodIx)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
v = zeros(d.G.cells.num, nsteps);
if isempty(d.currentDiagnostics(4).values)
    d = extractTracerForward(d, ts, injIx, prodIx);
end
tr = d.currentDiagnostics(4).values;
for k = 1:nsteps
    tr_cur  = tr(:, k);
    tof = sum(diagn(ts(k)).D.itracer(:, injIx).*diagn(ts(k)).D.itof(:, injIx),2);
    isfnt = tr_cur > 1e-5;
    tof(isfnt) = tof(isfnt)./tr_cur(isfnt);
    tof(~isfnt) = d.maxTOF;
    v(:, k) = min(tof, d.maxTOF);
end
d.currentDiagnostics(1).values = v/year;
end

function d = extractFAForward(d, ts, injIx, prodIx)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
v = zeros(d.G.cells.num, nsteps);
for k = 1:nsteps
    fa = min(diagn(ts(k)).D.ifa(:, injIx), [], 2);
    v(:, k) = min(fa, d.maxTOF);
end
d.currentDiagnostics(9).values = v/year;
end


function d = extractTOFBackward(d, ts, injIx, prodIx)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
v = zeros(d.G.cells.num, nsteps);
if isempty(d.currentDiagnostics(5).values)
    d = extractTracerBackward(d, ts, injIx, prodIx);
end
tr = d.currentDiagnostics(5).values;
for k = 1:nsteps
    tr_cur  = tr(:, k);
    tof = sum(diagn(ts(k)).D.ptracer(:, prodIx).*diagn(ts(k)).D.ptof(:, prodIx),2);
    isfnt = tr_cur > 1e-5;
    tof(isfnt) = tof(isfnt)./tr_cur(isfnt);
    tof(~isfnt) = d.maxTOF;
    v(:, k) = min(tof, d.maxTOF);
end
d.currentDiagnostics(2).values = v/year;
end

function d = extractFABackward(d, ts, injIx, prodIx)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
v = zeros(d.G.cells.num, nsteps);
for k = 1:nsteps
    fa = min(diagn(ts(k)).D.pfa(:, prodIx), [], 2);
    v(:, k) = min(fa, d.maxTOF);
end
d.currentDiagnostics(10).values = v/year;
end

function d = extractResidenceTime(d, ts, injIx, prodIx)
if isempty(d.currentDiagnostics(1).values)
    d = extractTOFForward(d, ts, injIx, prodIx);
end
if isempty(d.currentDiagnostics(2).values)
    d = extractTOFBackward(d, ts, injIx, prodIx);
end
trf = d.currentDiagnostics(1).values;
trb = d.currentDiagnostics(2).values;

d.currentDiagnostics(3).values = trf+trb;
end
    
function d = extractTracerProduct(d, ts, injIx, prodIx)
if isempty(d.currentDiagnostics(4).values)
    d = extractTracerForward(d, ts, injIx, prodIx);
end
if isempty(d.currentDiagnostics(5).values)
    d = extractTracerBackward(d, ts, injIx, prodIx);
end
trf = d.currentDiagnostics(4).values;
trb = d.currentDiagnostics(5).values;

d.currentDiagnostics(6).values = trf.*trb;
end

function d = extractSweepRegion(d, ts)
diagn  = d.Data.diagnostics;
nsteps = numel(ts);
p = 10; bc = 0.85;
cmap = d.Data.injColors;
% put last (reservoir-) color first since below we use D.ipart+1 (and the reservoir partition has value 0)
cmap = [cmap(end,:); cmap(1:end-1,:)];
rgb  = zeros(d.G.cells.num,3);
for k=1:nsteps
   maxconc = max(diagn(ts(k)).D.itracer, [], 2);
   w   = min(max(2 * (1 - maxconc), 0), 1);
   rgb = rgb + ...
      bsxfun(@plus, bsxfun(@times, 1-w.^p, cmap(diagn(ts(k)).D.ipart+1, :)), ...
                w.^p .* bc);
end
d.currentDiagnostics(7).values = rgb./nsteps;
end

function d = extractDrainageRegion(d, ts)
diagn = d.Data.diagnostics;
nsteps = numel(ts);
p = 10; bc = 0.85;
cmap = d.Data.prodColors;
% put reservoir-color first
cmap = [cmap(end,:); cmap(1:end-1,:)];
rgb  = zeros(d.G.cells.num,3);
for k=1:nsteps
   maxconc = max(diagn(ts(k)).D.ptracer, [], 2);
   w    = min(max(2 * (1 - maxconc), 0), 1);
   rgb = rgb + ...
      bsxfun(@plus, bsxfun(@times, 1-w.^p, cmap(diagn(ts(k)).D.ppart+1, :)), ...
                w.^p .* bc);
end
d.currentDiagnostics(8).values = rgb./nsteps;
end
