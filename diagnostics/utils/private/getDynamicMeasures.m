function m = getDynamicMeasures(d, tsel, wsel)
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

nT = numel(tsel.ix);
ninterp  = 51;
m.Phi = linspace(0,1,ninterp);
m.Ev  = m.Phi.^.25;
if isempty(wsel.injectorIx) && isempty(wsel.producerIx)
    % Compute Lorenz coefficient for the whole field for each
    % selected time step
    m.Ft  = zeros(ninterp,nT);
    m.tDt = zeros(ninterp,nT);
    m.LCt = zeros(nT,1);
    for i=1:nT
        D = d.Data.diagnostics(tsel.ix(i)).D;
        [F,Phi]   = computeFandPhi(d.G.cells.PORV,D.tof);
        [~, ii]   = unique(Phi);
        m.Ft(:,i) = interp1(Phi(ii), F(ii), m.Phi);
        m.LCt(i)  = computeLorenz(F, Phi);
        [Ev,tD]   = computeSweep(F, Phi);
        [~, ii]   = unique(Ev);
        m.tDt     = interp1(Ev(ii), tD(ii), m.Ev);
    end
    m.computePairs = false;
    m.wellName = [];
    return
    
elseif isempty(wsel.injectorIx)
    % Find all injectors that are communicating with the
    % selected producers
    prodIx = wsel.producerIx;
    injIx  = wsel.getConnectedInj(prodIx);
    
elseif isempty(wsel.producerIx)
    % Find all producers that are communicating with the selected
    % injectors
    injIx  = wsel.injectorIx;
    prodIx = wsel.getConnectedProd(injIx);
    
else
    injIx  = wsel.injectorIx;
    prodIx = wsel.producerIx;
end

[nI, nP] = deal(numel(injIx), numel(prodIx));
m.LCt = nan(nT,1);
m.Ft  = nan(ninterp,nT);
m.tDt = nan(ninterp,nT);
if (nI>1) && (nP>1)
    m.computePairs = false;
    m.wellName = [];
elseif nI==1
    m.LC  = nan(nP, nT);
    m.computePairs = true;
    m.names= arrayfun(@(x) x.label.String, ...
        d.WellPlot.producers(prodIx),'UniformOutput',false);
    m.wellName = d.WellPlot.injectors(injIx).label.String;
else
    m.LC  = nan(nI, nT);
    m.computePairs = true;
    m.names= arrayfun(@(x) x.label.String, ...
        d.WellPlot.injectors(injIx),'UniformOutput',false);
    m.wellName = d.WellPlot.producers(prodIx).label.String;
end

% Compute F-Phi and Lorenz coefficient for the whole region and
% for each well pair
[m.F, m.tD] = deal(nan(ninterp,nI*nP, nT));
for t=1:nT
    D = d.Data.diagnostics(tsel.ix(t)).D;
    
    [ritr, ritof, rptr, rptof] = deal(zeros(d.G.cells.num,1));
    n  = 1;
    for i=injIx
        itr   = D.itracer(:,i);
        itof  = D.itof(:,i);
        ritr  = ritr + itr;
        ritof = ritof + itr.*itof;
        for p=prodIx
            ptr   = D.ptracer(:,p);
            ptof  = D.ptof(:,p);
            rptr  = rptr + ptr;
            rptof = rptof + ptr.*ptof;
            if ~m.computePairs, continue, end
            ix = itr.*ptr>wsel.threshold;
            if sum(ix)>0
                [F,Phi] = computeFandPhi(d.G.cells.PORV(ix).*itr(ix).*ptr(ix), ...
                    [itof(ix) ptof(ix)]);
                [~, ii]     = unique(Phi);
                m.F(:,n,t)  = interp1(Phi(ii), F(ii), m.Phi);
                m.LC(n,t)   = computeLorenz(F, Phi);
                [Ev,tD]     = computeSweep(F, Phi);
                [~, ii]     = unique(Ev);
                m.tD(:,n,t) = interp1(Ev(ii), tD(ii), m.Ev);
            end
            n=n+1;
        end
    end
    
    ireg = ritr.*rptr;
    ix   = ireg > wsel.threshold;
    if nnz(ix) > 0
        [F,Phi] = computeFandPhi(d.G.cells.PORV(ix).*ireg(ix),...
            [ritof(ix)./ritr(ix) rptof(ix)./rptr(ix)]);
        [~, ii]    = unique(Phi);
        m.Ft(:,t)  = interp1(Phi(ii), F(ii), m.Phi);
        m.LCt(t)   = computeLorenz(F, Phi);
        [Ev,tD]    = computeSweep(F, Phi);
        [~, ii]    = unique(Ev);
        m.tDt(:,t) = interp1(Ev(ii), tD(ii), m.Ev);
    else
        m.LCt(t)  = nan;
    end
end
end
