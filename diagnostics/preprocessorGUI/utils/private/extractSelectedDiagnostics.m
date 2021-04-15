function [d, vals, lims, flag] = extractSelectedDiagnostics(d, prop, modsel, wsel)
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

ms = modsel.ix;
if isempty(ms)
    ms = 1;
end


[injIx, prodIx] = deal(wsel.injectorIx, wsel.producerIx);
if isempty(injIx)
   injIx = 1:numel(d.Data{ms(1)}.diagnostics.D.inj);
end
if isempty(prodIx)
   prodIx = 1:numel(d.Data{ms(1)}.diagnostics.D.prod);
end
flag = false; % true if 'vals' is a color a not a set of scalar values

switch prop
    case d.currentDiagnostics(1).name % TOF forward
        if isempty(d.currentDiagnostics(1).values)
            d = extractTOFForward(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(1).values;
        lims = [0 min(max(vals(:)), d.maxTOF/year)];
        
    case d.currentDiagnostics(2).name % TOF backward
        if isempty(d.currentDiagnostics(2).values)
            d = extractTOFBackward(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(2).values;
        lims = [0 min(max(vals(:)), d.maxTOF/year)];
        
    case d.currentDiagnostics(3).name % Residence time
        if isempty(d.currentDiagnostics(3).values)
            d = extractResidenceTime(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(3).values;
        lims = [0 min(max(vals(:)),2*d.maxTOF/year)];
        
    case d.currentDiagnostics(4).name % Tracer forward
        if isempty(d.currentDiagnostics(4).values)
            d = extractTracerForward(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(4).values;
        lims = [0 1];
        
    case d.currentDiagnostics(5).name % Tracer backward
        if isempty(d.currentDiagnostics(5).values)
            d = extractTracerBackward(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(5).values;
        lims = [0 1];
        
    case d.currentDiagnostics(6).name % Tracer product
        if isempty(d.currentDiagnostics(6).values)
            d = extractTracerProduct(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(6).values;
        lims = [0 1];

    case d.currentDiagnostics(7).name % Sweep region
        if isempty(d.currentDiagnostics(7).values)
           d = extractSweepRegion(d, ms);
        end
        vals = d.currentDiagnostics(7).values;
        lims = [0 1];
        flag = true;

    case d.currentDiagnostics(8).name % Drainage region
        if isempty(d.currentDiagnostics(8).values)
           d = extractDrainageRegion(d, ms);
        end
        vals = d.currentDiagnostics(8).values;
        lims = [0 1];
        flag = true;
        
    case d.currentDiagnostics(9).name % First arrival forward
        if isempty(d.currentDiagnostics(9).values)
           d = extractFAForward(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(9).values;
        lims = [0 d.maxTOF/year];
        
    case d.currentDiagnostics(10).name % First arrival backward
        if isempty(d.currentDiagnostics(10).values)
           d = extractFABackward(d, ms, injIx, prodIx);
        end
        vals = d.currentDiagnostics(10).values;
        lims = [0 d.maxTOF/year];
        
    otherwise
        warning(['Diagnostics property ', prop, ' is not known or badly spelled....'])
end
end

function d = extractTracerForward(d, ms, injIx, ~)
nmods = numel(ms);
v = zeros(d.Data{1}.G.cells.num, nmods);
for k = 1:nmods
    diagn = d.Data{ms(k)}.diagnostics;
    v(:, k) = sum(diagn.D.itracer(:, injIx), 2);
end
d.currentDiagnostics(4).values = v;
end

function d = extractTracerBackward(d, ms,  ~, prodIx)
nmods = numel(ms);
v = zeros(d.Data{1}.G.cells.num, nmods);
for k = 1:nmods
    diagn = d.Data{ms(k)}.diagnostics;
    v(:, k) = sum(diagn.D.ptracer(:, prodIx), 2);
end
d.currentDiagnostics(5).values = v;
end

function d = extractTOFForward(d, ms, injIx, prodIx)

nmods = numel(ms);
v = zeros(d.Data{1}.G.cells.num, nmods);
if isempty(d.currentDiagnostics(4).values)
    d = extractTracerForward(d, ms, injIx, prodIx);
end
tr = d.currentDiagnostics(4).values;
for k = 1:nmods
    diagn = d.Data{ms(k)}.diagnostics;
    tr_cur  = tr(:, k);
    tof = sum(diagn.D.itracer(:, injIx).*diagn.D.itof(:, injIx),2);
    isfnt = tr_cur > 1e-5;
    tof(isfnt) = tof(isfnt)./tr_cur(isfnt);
    tof(~isfnt) = d.maxTOF;
    v(:, k) = min(tof, d.maxTOF);
end
d.currentDiagnostics(1).values = v/year;
end

function d = extractFAForward(d, ms, injIx, prodIx)
nmods = numel(ms);
v = zeros(d.Data{1}.G.cells.num, nmods);
for k = 1:nmods
    diagn = d.Data{ms(k)}.diagnostics;    
    fa = min(diagn.D.ifa(:, injIx), [], 2);
    v(:, k) = min(fa, d.maxTOF);
end
d.currentDiagnostics(9).values = v/year;
end


function d = extractTOFBackward(d, ms, injIx, prodIx)
nmods = numel(ms);
v = zeros(d.Data{1}.G.cells.num, nmods);
if isempty(d.currentDiagnostics(5).values)
    d = extractTracerBackward(d, ms, injIx, prodIx);
end
tr = d.currentDiagnostics(5).values;
for k = 1:nmods
    diagn = d.Data{ms(k)}.diagnostics;
    tr_cur  = tr(:, k);
    tof = sum(diagn.D.ptracer(:, prodIx).*diagn.D.ptof(:, prodIx),2);
    isfnt = tr_cur > 1e-5;
    tof(isfnt) = tof(isfnt)./tr_cur(isfnt);
    tof(~isfnt) = d.maxTOF;
    v(:, k) = min(tof, d.maxTOF);
end
d.currentDiagnostics(2).values = v/year;
end


function d = extractFABackward(d, ms, injIx, prodIx)
nmods = numel(ms);
v = zeros(d.Data{1}.G.cells.num, nmods);
for k = 1:nmods
    diagn = d.Data{ms(k)}.diagnostics;    
    fa = min(diagn.D.pfa(:, prodIx), [], 2);
    v(:, k) = min(fa, d.maxTOF);
end
d.currentDiagnostics(10).values = v/year;
end


function d = extractResidenceTime(d, ms, injIx, prodIx)
if isempty(d.currentDiagnostics(1).values)
    d = extractTOFForward(d, ms, injIx, prodIx);
end
if isempty(d.currentDiagnostics(2).values)
    d = extractTOFBackward(d, ms, injIx, prodIx);
end
trf = d.currentDiagnostics(1).values;
trb = d.currentDiagnostics(2).values;

d.currentDiagnostics(3).values = trf+trb;
end
    
function d = extractTracerProduct(d, ms, injIx, prodIx)
if isempty(d.currentDiagnostics(4).values)
    d = extractTracerForward(d, ms, injIx, prodIx);
end
if isempty(d.currentDiagnostics(5).values)
    d = extractTracerBackward(d, ms, injIx, prodIx);
end
trf = d.currentDiagnostics(4).values;
trb = d.currentDiagnostics(5).values;

d.currentDiagnostics(6).values = trf.*trb;
end

function d = extractSweepRegion(d, ms)
nmods = numel(ms);
rgb  = zeros(d.Data{1}.G.cells.num,3);    
for k=1:nmods
    diagn = d.Data{ms(k)}.diagnostics;
    p = 10; bc = 0.85;
    cmap = d.Data{ms(k)}.injColors;
    % put last (reservoir-) color first since below we use D.ipart+1 (and the reservoir partition has value 0)
    cmap = [cmap(end,:); cmap(1:end-1,:)];
    maxconc = max(diagn.D.itracer, [], 2);
    w   = min(max(2 * (1 - maxconc), 0), 1);
    rgb = rgb + ...
      bsxfun(@plus, bsxfun(@times, 1-w.^p, cmap(diagn.D.ipart+1, :)), ...
                w.^p .* bc);
end
d.currentDiagnostics(7).values = rgb./nmods;
end

function d = extractDrainageRegion(d, ms)

nmods = numel(ms);
rgb  = zeros(d.Data{1}.G.cells.num,3);
for k=1:nmods
    diagn = d.Data{ms(k)}.diagnostics;
    p = 10; bc = 0.85;
    cmap = d.Data{ms(k)}.prodColors;
    % put reservoir-color first
    cmap = [cmap(end,:); cmap(1:end-1,:)];

    maxconc = max(diagn.D.ptracer, [], 2);
    w    = min(max(2 * (1 - maxconc), 0), 1);
    rgb = rgb + ...
        bsxfun(@plus, bsxfun(@times, 1-w.^p, cmap(diagn.D.ppart+1, :)), ...
        w.^p .* bc);
end
d.currentDiagnostics(8).values = rgb./nmods;
end
