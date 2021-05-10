function SCHEDULE = convertScheduleToDeck(model, schedule, varargin)
% Create deck-type schedule-structure from MRST model and schedule. 
%
% SYNOPSIS:
%   SCHEDULE = convertScheduleToDeck(model, schedule)
%   SCHEDULE = convertScheduleToDeck(model, schedule, 'pn1', pv1, ...)
%
% DESCRIPTION:
%   This function creates a deck-type schedule-structure amenable for 
%   writing input files to other simulators by e.g., writeSchedule.
%
% REQUIRED PARAMETERS:
%   model     -  MRST model  
%
%   schedule  - MRST schedule
%
% OPTIONAL PARAMETERS:
%   'linearIndex'   - If true, assumes grid is unstructured and output as 1D, 
%                     hence cell n coresponds to [I,J,K] = [n, 1, 1]. If
%                     false, [I,J,K] is deduced from G.cartdims. 
%                     Default value: true.
%   'reduceOutput'  - If true, write COMPDAT, WELSPECS, WCONINJE and WCONPROD 
%                     for a given control step only if changed from previous 
%                     control step. If false, all data are writted for each 
%                     step. Default value: true.

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

opt = struct('linearIndex', true, ...
             'reduceOutput', true);
opt = merge_options(opt,varargin{:});
% prototypes
% COMPDAT       1        2  3  4  5  6       7    8  9       10  11  12        13   14
%               nm       I  J  K1 K2 flag  satnum wi d       Kh  skin D-fac    dir  r0       
COMPDAT_tmp = {'INJE01', 1, 1, 1, 1, 'OPEN', -1, -1, 0.1000, -1, 0, 'Default', 'Z', -1;...
               'PROD01', 3, 1, 1, 1, 'OPEN', -1, -1, 0.1000, -1, 0, 'Default', 'Z', -1};

% WELSPECS        1       2    3  4  5     6       7  8      9       10    11 12     13
%                nm      grp  I  J  dep   pph  drain infl   shutin  cr-fl pvtn dens  fip                                                                       
WELSPECS_tmp = {'INJE01', 'I', 1, 1, NaN, 'WATER', 0, 'STD', 'SHUT', 'YES', 0, 'SEG', 0;
                'PROD01', 'P', 1, 1, NaN ,'OIL',   0 ,'STD', 'SHUT', 'YES', 0 ,'SEG', 0 };

v = Inf;
% WCONINJE      1         2        3       4       5    6    7   8   9   10  11 12 13 14    
%               nm        type     flag    cntr    rate resv bhp thp vfp rv  th ro rw rg               
WCONINJE_tmp = {'INJE01', 'WATER', 'SHUT', 'RATE', v,   v,   v,   v, 0,  v,  v, v, v, v};


% WCONPROD       1        2       3      4    5    6    7    8    9   10  11  12          
%                nm       flag    cntr   orat wrat grat lrat resv bhp thp vfp lift
WCONPROD_tmp = {'PROD01', 'SHUT', 'BHP', v,   v,   v,   v,   v,   v,  v,  0   v};

%
% WCONHIST       1        2       3      4    5    6    7    8    9   10  11  12          
%                nm       flag    cntr   orat wrat grat vfp vlfq thp bhp  
WCONHIST_tmp = {'PROD01', 'OPEN', 'BHP', v,   v,   v,   0  , 0,  v, ...
                v};

%nonre
%GCONINJE

%
% phases
phaseNames = num2cell(model.getPhaseNames);
for k = 1:numel(phaseNames)
    switch phaseNames{k}
        case 'W'
            phaseNames{k} = 'WATER';
        case 'O'
            phaseNames{k} = 'OIL';
        case 'G'
            phaseNames{k} = 'GAS';
    end
end

reduce = opt.reduceOutput;

SCHEDULE.RPTSCHED = 1;
% steps                   
SCHEDULE.step = struct('control', schedule.step.control, ...
                       'val',     schedule.step.val);                   
% IJK-map
nc = model.G.cells.num;
if opt.linearIndex
    cartDims = [nc, 1, 1];
    indexMap = (1:nc)';
else % map model.G to box
    cartDims = model.G.cartDims;
    indexMap = model.G.cells.indexMap;
end
getIJK = @(c)ind2sub(cartDims, indexMap(c));

nSteps  = numel(schedule.control);
cntr    = struct('WELSPECS', {{}}, 'COMPDAT', {{}}, 'WCONINJE', {{}}, 'WCONPROD', {{}});
control = repmat(cntr, nSteps, 1);

for sno = 1:numel(control)
    % WELSPECS ------------------------------------------------------------
    W = schedule.control(sno).W;
    ix = getUpdateIx(schedule, sno, 'WELSPECS', reduce);
    WELSPECS = cell(numel(ix), size(WELSPECS_tmp, 2));
    for ixno = 1:numel(ix)
        wno = ix(ixno);
        if W(wno).sign > 0 %inj
            WELSPECS(ixno,:) = WELSPECS_tmp(1,:);
            [~, pphix] = max(W(wno).compi);
            WELSPECS{ixno,6} = phaseNames{pphix};
        else
            WELSPECS(ixno,:) = WELSPECS_tmp(2,:);
        end
        WELSPECS{ixno,1} = W(wno).name;
        [I, J, ~] = getIJK(W(wno).cells(1));
        WELSPECS{ixno,3} = I;
        WELSPECS{ixno,4} = J;
        WELSPECS{ixno,5} = W(wno).refDepth;
        % preferred phase
    end
    control(sno).WELSPECS = WELSPECS;

    % COMPDAT -------------------------------------------------------------
    ix = getUpdateIx(schedule, sno, 'COMPDAT', reduce);
    COMPDAT = cell(numel(W), 1);
    for wno = 1:numel(W)
        cix = ix{wno};
        cells = W(wno).cells(cix);
        if W(wno).sign > 0 %inj
            tmpl = COMPDAT_tmp(1,:);
        else
            tmpl = COMPDAT_tmp(2,:);
        end
        tmpl{1} = W(wno).name;
        COMPDAT{wno} = repmat(tmpl, [numel(cells), 1]);
        if ~isempty(cix)
            [I, J, K] = getIJK(cells);
            % ijk
            COMPDAT{wno}(:, 2:5) = num2cell([I, J, K, K]);
            % open/shut
            flags = W(wno).cstatus(cix);
            nms   = {'SHUT'; 'OPEN'};
            COMPDAT{wno}(:, 6) = nms(flags+1);
            % WI
            COMPDAT{wno}(:, 8) = num2cell( W(wno).WI(cix) );
            % r
            r = W(wno).r;
            if numel(r) == 1
                r = repmat(r(1), [numel(cells), 1]);
            else
                r = r(cix);
            end
            COMPDAT{wno}(:, 9) = num2cell(r*2);
            % dir
            COMPDAT{wno}(:, 13) = num2cell(upper(W(wno).dir(cix)));
        end
    end
    control(sno).COMPDAT = vertcat(COMPDAT{:});

    % WCONINJE ------------------------------------------------------------
    if isempty(W)
        injIx = [];
    else
        injIx = find([W.sign] > 0);
    end
    WCONINJE = repmat(WCONINJE_tmp, [numel(injIx), 1]);
    for ino = 1:numel(injIx)
        wno = injIx(ino);
        WCONINJE{ino,1} = W(wno).name;
        % deal with injection comp
        wog = [0 0 0];
        ph  = {'W', 'O', 'G'};
        for kp = 1:3
            ix = model.getPhaseIndex(ph{kp});
            if ~isempty(ix)
                wog(kp) = W(wno).compi(ix);
            end
        end
        itypes = {'WATER', 'OIL', 'GAS'};
        WCONINJE{ino,2} = itypes{find(wog, 1, 'first')};
        if nnz(wog) > 1
            WCONINJE(ino, 12:14) = num2cell(wog([2 1 3]));
        end
        % open/shut
        if(W(wno).status)
            WCONINJE{ino,3} = 'OPEN';
        else
            WCONINJE{ino,3} = 'SHUT';
        end
        % control mode and limits
        WCONINJE{ino,4} = upper(W(wno).type);
        lims = W(wno).lims;
        % possibly overwrite "target limit"
        lims.(W(wno).type) = W(wno).val;
        fnms = fieldnames(lims);
        for kn = 1:numel(fnms)
            fnm = fnms{kn};
            if isfinite(lims.(fnm))
                switch fnm
                    case 'bhp'
                        WCONINJE{ino,7} = lims.(fnm);
                    case 'thp'
                        WCONINJE{ino,8} = lims.(fnm);
                    case 'resv'
                        WCONINJE{ino,6} = abs(lims.(fnm));
                    case 'rate'
                        WCONINJE{ino,5} = abs(lims.(fnm));
                    otherwise
                        warning([fnm, ' limit present for injector '])
                        %error('Not done');
                end
            end
        end
        %end
    end
    control(sno).WCONINJE = WCONINJE;

    % WCONPROD ------------------------------------------------------------
    if isempty(W)
        prodIx = [];
    else
        prodIx = [W.sign] <= 0;
    end
    prodIx_pred=prodIx;
    prodIx_hist=prodIx;
    for pno=1:numel(prodIx)
        if(prodIx(pno))
        wno = pno;%prodIx(pno);
        is_history= strcmp(lower(W(wno).type),'resv_history');
        if(is_history)
            prodIx_pred(wno)=false;
        else
            prodIx_hist(wno)=false;
        end
        end
        if(not(W(pno).status))
            prodIx_pred(wno)=false;
            prodIx_hist(wno)=false;
        end
    end
    prodIx_pred=find(prodIx_pred);
    prodIx_hist=find(prodIx_hist);
    WCONPROD = makeWCONPROD(WCONPROD_tmp,prodIx_pred,W);
    WCONHIST = makeWCONHIST(WCONHIST_tmp,prodIx_hist,W);
    control(sno).WCONPROD = WCONPROD;
    control(sno).WCONHIST = WCONHIST;
end

SCHEDULE.control = control;
end

%--------------------------------------------------------------------------

function WCONPROD = makeWCONPROD(WCONPROD_tmp,prodIx,W)
        
    WCONPROD = repmat(WCONPROD_tmp, [numel(prodIx), 1]);
    pno=0;
    for pno = 1:numel(prodIx)
        wno = prodIx(pno);
        WCONPROD{pno,1} = W(wno).name;
        % open/shut
        if W(wno).status
            WCONPROD{pno,2}='OPEN';
        else
            WCONPROD{pno,2}='SHUT';
        end
        % control mode and limits
        WCONPROD{pno,3} = upper(W(wno).type);
        lims = W(wno).lims;
        if isempty(lims)
            continue
        end
        % possibly overwrite "target limit"
        lims.(W(wno).type) = W(wno).val;
        fnms = fieldnames(lims);
        for kn = 1:numel(fnms)
            fnm = fnms{kn};
            if isfinite(lims.(fnm))
                switch lower(fnm)
                    case 'bhp'
                        WCONPROD{pno,9} = lims.(fnm);
                    case 'thp'
                        WCONPROD{pno,10} = lims.(fnm);
                    case {'resv'}
                        WCONPROD{pno,8} = abs(lims.(fnm));
                    case 'orat'
                        WCONPROD{pno,4} = abs(lims.(fnm));
                    case 'grat'
                        WCONPROD{pno,6} = abs(lims.(fnm));
                    case 'wrat'
                        WCONPROD{pno,5} = abs(lims.(fnm));
                    case 'lrat'
                        WCONPROD{pno,7} = abs(lims.(fnm));
                    case 'grup'
                    % this may happen if inactive wells is not removed
                    otherwise
                        error('Not done');
                end
            end
        end
    end
end

%--------------------------------------------------------------------------

function WCONHIST = makeWCONHIST(WCONHIST_tmp,prodIx,W)
        
    WCONHIST = repmat(WCONHIST_tmp, [numel(prodIx), 1]);
    pno=0;
    for pno = 1:numel(prodIx)
        wno = prodIx(pno);
        WCONHIST{pno,1} = W(wno).name;
        % open/shut
        if W(wno).status
            WCONHIST{pno,2}='OPEN';
        else
            WCONHIST{pno,2}='SHUT';
        end
        % control mode and limits
        type =  upper(W(wno).type);
        if(strcmp(lower(type),'resv_history'))
            type = 'resv';
            %elseif(type == 'lrat_history')
            %    type = 'lrat';
            %elseif(type == 'orat_history')
            %    type = 'orat';   
        else
            error(['This is not a histry mode', type] );
        end
        vals = W(wno).compi*W(wno).val;
        WCONHIST{pno,3} = upper(type);
        lims = W(wno).lims;
        lims.grat=-vals(3);
        lims.wrat=-vals(1);
        lims.orat=-vals(2);
        fnms = fieldnames(lims);
        for kn = 1:numel(fnms)
            fnm = fnms{kn};
            if isfinite(lims.(fnm))
                switch lower(fnm)
                    case 'bhp'
                        WCONHIST{pno,10} = lims.(fnm);
                    case 'thp'
                        WCONHIST{pno,9} = lims.(fnm);
                    case 'orat'
                        WCONHIST{pno,4} = abs(lims.(fnm));
                    case 'grat'
                        WCONHIST{pno,6} = abs(lims.(fnm));
                    case 'wrat'
                        WCONHIST{pno,5} = abs(lims.(fnm));
                    otherwise
                        error('Not done');
                end
            end
        end
    end
end

%--------------------------------------------------------------------------

function ix = getUpdateIx(schedule, step, fld, doReduce)
W  = schedule.control(step).W;
[nw, nconn] = deal(numel(W), arrayfun(@(w)numel(w.cells), W));
if step > 1
    Wp = schedule.control(step-1).W;
end
switch fld
    case 'WELSPECS'
        ix = (1:nw)';
        if step > 1 && doReduce
            ix = find(vertcat(W.sign)~=vertcat(Wp.sign));
        end
    case 'COMPDAT'
        ix = cellfun(@(n)(1:n)', num2cell(nconn), 'UniformOutput', false);
        if step > 1 && doReduce
            for kw = 1:nw
                ix{kw} = find( (W(kw).cstatus ~= Wp(kw).cstatus) | (W(kw).WI ~= Wp(kw).WI) );         
            end
        end
end
end
