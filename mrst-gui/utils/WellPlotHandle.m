classdef WellPlotHandle < handle
    properties
        injectors
        producers
        injectorNames
        producerNames
    end
    properties (SetAccess = immutable)
        nInj
        nProd
        nCases 
        map
    end
    properties (Dependent, SetObservable)
        Visible              % 'on' / 'off'
        visibleInjectors     % logic array of size nInj x 1  
        visibleProducers     % logic array of size nProd x 1  
        visibleCases         % logic array of size nCases x 1  
        injectorColor
        producerColor
        closedColor
    end
    % Convenience props to aid functionality of set/get props
    properties (Hidden)
        caseIx  = 1;
        injIx
        prodIx
        iColor  = [.2 .2 .7];
        pColor  = [.7 .2 .2];
        cColor  = [.7 .7 .7];
        statusExtended
        Parent
    end
    
    methods
        
        function s = WellPlotHandle(G, W, varargin)
            opt = struct('Parent',     [], ...
                         'Visible',                      'on', ...      % all wells visible
                         'visibleCases',                     1, ...
                         'injectorColor',           [.2 .2 .7], ...
                         'producerColor',           [.7 .2 .2], ...
                         'closedColor',             [.7 .7 .7], ...     % default grey
                         'distinguishConnectionStatus',  false); 
            [opt, extraOpt] = merge_options(opt, varargin{:});
              
            if isempty(opt.Parent)
                opt.Parent = gca;
            end
            
            % treat as multiple well-cases
            if ~iscell(W)
                W = {W};
            end
            s.Parent = opt.Parent;
            s.nCases = numel(W);
            isInj    = [W{1}.sign]>0;
            s.nInj   = nnz(isInj);
            s.nProd  = nnz(~isInj);
            s.injectorNames = arrayfun(@(x)x.name, W{1}(isInj), 'UniformOutput', false);
            s.producerNames = arrayfun(@(x)x.name, W{1}(~isInj), 'UniformOutput', false);
                        
            [Ws, s.map] = getWellSuperset(W, opt.distinguishConnectionStatus);

            if ~isfield(Ws(1), 'trajectory')
                Ws = addTrajectories(Ws, G, 10, opt.distinguishConnectionStatus);
            end
            
            isInj = [Ws.sign]>0;
            s.injectors  = plotWellData(G, Ws(isInj), extraOpt{:}, ...
                'linePlot', true, 'color', s.iColor, ...
                'labelBackgroundColor',[.75 .75 .75],'labelFontSize',8);
            s.injectors  = setParents(s.injectors(:), opt.Parent);
            set([s.injectors.label], 'interpreter', 'none')
            s.statusExtended.inj = true(numel(s.injectors),1);
            
            
            s.producers = plotWellData(G, Ws(~isInj), extraOpt{:}, ...
                'linePlot', true, 'color', s.pColor, ...
                'labelBackgroundColor',[.75 .75 .75],'labelFontSize',8);
            s.producers = setParents(s.producers(:), opt.Parent);
            set([s.producers.label], 'interpreter', 'none')
            s.statusExtended.prod = true(numel(s.producers),1);
            
            % set initial colors
            s.iColor        = opt.injectorColor;
            s.pColor        = opt.producerColor;
            s.cColor        = opt.closedColor;   
            s.Visible       = opt.Visible; 
            s.visibleCases  = opt.visibleCases;
        end
        %------------------------------------------------------------------
        function set.injectorColor(s, val)
            fn = fieldnames(s.injectors);
            for k = 1:numel(fn)
                set([s.injectors.(fn{k})], 'Color', val);
            end
            s.iColor = val;
            s.closedColor = s.closedColor;
        end
        function val = get.injectorColor(s)
            val = s.iColor;
        end
        %------------------------------------------------------------------
        function set.producerColor(s, val)
            fn = fieldnames(s.producers);
            for k = 1:numel(fn)
                set([s.producers.(fn{k})], 'Color', val);
            end
            s.pColor = val;
            s.closedColor = s.closedColor;
        end
        function val = get.producerColor(s)
            val = s.pColor;
        end
        %------------------------------------------------------------------
        function set.closedColor(s, val)
            s.cColor = val;
            s.statusExtended.inj = true(size(s.injectors));
            s.visibleInjectors = s.visibleInjectors;
            s.statusExtended.prod = true(size(s.producers));
            s.visibleProducers    = s.visibleProducers;
        end
        function val = get.closedColor(s)
            val = s.cColor;
        end
        %------------------------------------------------------------------
        function set.Visible(s, val)
            s.visibleInjectors = val;
            s.visibleProducers = val;
        end
        function val = get.Visible(s)
            if all([s.visibleInjectors; s.visibleProducers])
                val = 'on';
            elseif all([~s.visibleInjectors; ~s.visibleProducers])
                val = 'off';
            else
                val = 'mixed';
            end
        end
        %------------------------------------------------------------------
        function set.visibleInjectors(s, val)
            if ischar(val)
                if strcmp(val, 'on')
                    val = true(s.nInj, 1);
                else
                    val = false(s.nInj, 1);
                end
            elseif ~islogical(val)
                tmp = false(s.nInj, 1);
                tmp(val) = true;
                val = tmp;
            end
            if numel(val) ~= s.nInj
                error('Index exceeds the number of injectors (%d)', s.nInj);
            end
            s.injIx = val;
            [val, stat] = s.sub2super(val, 'inj');
            ix = val~=s.getVisibleInjectorsExtended();
            if any(ix)
                ix = find(ix);      
                fn = fieldnames(s.injectors);
                opts = {'off', 'on'};
                for k = 1:numel(ix)
                    for n = 1:numel(fn)
                        set([s.injectors(ix(k)).(fn{n})], 'Visible', opts{val(ix(k))+1})
                    end
                end
            end
            s.updateClosedInjectors(stat);
        end
        function val = get.visibleInjectors(s)
            val = s.injIx;
        end
        %------------------------------------------------------------------    
        function set.visibleProducers(s, val)
            if ischar(val)
                if strcmp(val, 'on')
                    val = true(s.nProd, 1);
                else
                    val = false(s.nProd, 1);
                end
            elseif ~islogical(val)
                tmp = false(s.nProd, 1);
                tmp(val) = true;
                val = tmp;
            end
            if numel(val) ~= s.nProd
                error('Index exceeds the number of producers (%d)', s.nProd);
            end
            s.prodIx = val;
            [val, stat] = s.sub2super(val, 'prod');
            ix = val~=s.getVisibleProducersExtended();
            if any(ix)
                ix = find(ix);            
                fn = fieldnames(s.producers);
                opts = {'off', 'on'};
                for k = 1:numel(ix)
                    for n = 1:numel(fn)
                        set([s.producers(ix(k)).(fn{n})], 'Visible', opts{val(ix(k))+1})
                    end
                end
            end
            s.updateClosedProducers(stat);
        end
        function val = get.visibleProducers(s)
            val = s.prodIx;
        end
        %------------------------------------------------------------------
        function set.visibleCases(s, val)
            if ischar(val)
                if strcmp(val, 'on')
                    val = true(s.nCases, 1);
                else
                    val = false(s.nCases, 1);
                end
            elseif ~islogical(val)
                tmp = false(s.nCases, 1);
                tmp(val) = true;
                val = tmp;
            end
            if numel(val) ~= s.nCases
                error('Index exceeds the number of well-cases (%d)', s.nCases);
            end
            s.caseIx = val;
            % reset
            s.visibleInjectors = s.injIx;
            s.visibleProducers = s.prodIx;
        end
        function val = get.visibleCases(s)
            val = s.caseIx;
        end

        %------------------------------------------------------------------
        
        function [ixs, stat] = sub2super(s, ix, tp)
            if strcmp(tp, 'inj')
                ixs = false(numel(s.injectors),1);
            else
                ixs = false(numel(s.producers),1);
            end
            k = setdiff( unique(s.map.(tp)(ix, s.caseIx)), 0);
            ixs(k) = true;
            if nargout > 1
                k2   = setdiff( unique(s.map.open.(tp)(ix, s.caseIx)), 0);
                stat = false(size(ixs));
                stat(k2) = true;
            end
                
        end

        function ix = super2sub(s, ix, tp)
             if ~isempty(s.map)
                 m   = max(s.map.(tp), [], 2);
                 tmp = m(ix);
                 if strcmp(tp, 'inj')
                     ix  = false(s.nInj, 1);
                 else
                     ix = false(s.nProd, 1);
                 end
                 ix(tmp) = true;
            end
        end  
        
        function [] = updateClosedInjectors(s, stat)
            upd     =  stat ~= s.statusExtended.inj;
            upd     = upd & getVisibleInjectorsExtended(s);
            if any(upd)
                fn = {'connector', 'body', 'label'};
                for k = 1:numel(fn)
                    set([s.injectors(upd & ~stat).(fn{k})], 'Color', s.closedColor);
                    set([s.injectors(upd &  stat).(fn{k})], 'Color', s.injectorColor);
                end
                s.statusExtended.inj(upd) = ~s.statusExtended.inj(upd);
            end
        end
        
        function [] = updateClosedProducers(s, stat)
            upd     =  stat ~= s.statusExtended.prod;
            upd     = upd & getVisibleProducersExtended(s);
            if any(upd)
                fn = {'connector', 'body'};
                for k = 1:numel(fn)
                    set([s.producers(upd & ~stat).(fn{k})], 'Color', s.closedColor);
                    set([s.producers(upd &  stat).(fn{k})], 'Color', s.producerColor);
                end
                s.statusExtended.prod(upd) = ~s.statusExtended.prod(upd);
            end
        end
        %------------------------------------------------------------------
        function ix = getVisibleInjectorsExtended(s)
            ix = arrayfun(@(x)strcmp(x.label.Visible, 'on'), s.injectors);
        end
        
        function ix = getVisibleProducersExtended(s)
            ix = arrayfun(@(x)strcmp(x.label.Visible, 'on'), s.producers);
        end
        
        function delete(s)
            tp   = {'injectors', 'producers'};
            flds = {'label', 'connector', 'body'};
            for kt = 1:numel(tp)
                for kw = 1:numel(s.(tp{kt}))
                    for kf = 1:numel(flds)
                        delete(s.(tp{kt})(kw).(flds{kf}));
                    end
                end
            end
        end
    end
end

function wp = setParents(wp, p)
fn = fieldnames(wp);
for k = 1:numel(fn)
    set([wp.(fn{k})], 'Parent', p)
end
end
    
function W = addTrajectories(W, G, np, connOpt)
% add well trajectory for plotting, currently just np points along a single
% quadratic curve

for k = 1:numel(W)
    cst = W(k).cstatus;
    if ~connOpt
        cst = true(size(cst));
    end
    c  = G.cells.centroids(W(k).cells(cst), :);
    if isempty(c)
        c  = G.cells.centroids(W(k).cells(1), :);
    end
    dim = G.griddim;
    if dim < 3
        c = [c, zeros(size(c, 1), 3 - dim)];
    end
    if ~strcmp(W(k).type, 'aquifer')
        nc = size(c,1);
        c0 = c(1,:);
        % take point furthest away as toe
        v  = c-ones(nc,1)*c0;
        d2 = sum(v.*v, 2);
        [~, ix] = max(d2);
        c2 = c(ix,:);
        
        s = linspace(0, 1, np)';
        [f0, f1, f2] = deal((1-s).^2, 2*(1-s).*s, s.^2);
        % x,y,z - coords take second point as mean shifted twice away from c0-cn
        mc = mean(c);
        %
        n   = (c2-c0)/norm(c2-c0);
        pmc = c0 - dot(c0-mc, n)*n;
        c1  =  2*mc - pmc;
        c1(3) = min(c2(3),c1(3));
        traj  = f0*c0 + f1*c1 + f2*c2;
        
        
        %s = (linspace(0, 1, np).^(1/3))';
        %[f0, f1, f2] = deal((1-s).^2, 2*(1-s).*s, s.^2);
        %traj(:, 3) = f0*c0(3) + f1*c1(3) + f2*c2(3);
        
        W(k).trajectory = traj;
    else
        W(k).trajectory = c(1,:);
    end
end
end



function [W, map] = getWellSuperset(cases, distinguishConnectionStatus)
% well-names are assumed to be the same for each case (might want to check this)
nw    = numel(cases{1});
isInj = [cases{1}.sign] > 0;
[injNo, prodNo] = deal(cumsum(isInj), cumsum(~isInj));
[ni, np] = deal(nnz(isInj), nnz(~isInj)); 
nc = numel(cases);
uniqueWells    = cell(nw ,1);
uniqueNo       = cell(nw ,1);
[mapi, mapp]   = deal(nan(ni, nc), nan(np, nc));
status         = false(nw, nc);
[ci, cp]  = deal(0);
for kw = 1:nw
    for kc = 1:nc
        w = cases{kc}(kw);
        [uniqueWells{kw}, hit] = mergeWells(uniqueWells{kw}, w, distinguishConnectionStatus);
        if hit == 0
            if isInj(kw)
                ci = ci + 1;
                mapi(injNo(kw),kc) = ci;
                uniqueNo{kw} = [uniqueNo{kw}, ci];
            else
                cp = cp + 1;
                mapp(prodNo(kw),kc) = cp;
                uniqueNo{kw} = [uniqueNo{kw}, cp]; 
            end
        else
            if isInj(kw)
                mapi(injNo(kw),kc) = uniqueNo{kw}(hit);
            else
                mapp(prodNo(kw),kc) = uniqueNo{kw}(hit);
            end
        end
        status(kw,kc) = w.status && any(w.cstatus);
    end
end
% add suffix to mutiple well instances
for kw = 1:nw
    act = find( [uniqueWells{kw}.status] );
    if numel(act) > 1
        for ku = 1:numel(act)
            uniqueWells{kw}(act(ku)).name = sprintf('%s_[%s]', uniqueWells{kw}(act(ku)).name, char('A'+ku-1));
        end
    end
end

W = vertcat(uniqueWells{:});
map = struct('inj',       mapi,             'prod',      mapp, ...
             'open', struct('inj',  status(isInj,:).*mapi, ...
                            'prod', status(~isInj,:).*mapp));
end

function [W, hit] = mergeWells(W, w, connFlag)
hit = 0;
for k = 1:numel(W)
    if isSameWell(W(k), w, connFlag)
        hit = k;
        break;
    end
end
if hit == 0
    W = [W;w];
end
end

function b = isSameWell(w1, w2, connFlag)
b = false;
if numel(w1.cells) == numel(w2.cells)
    if connFlag
        b = all(w1.cstatus.*w1.cells == w2.cstatus.*w2.cells);
    else
        b = all(w1.cells == w2.cells);
    end
end
end

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
