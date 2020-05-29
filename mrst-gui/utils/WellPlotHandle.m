classdef WellPlotHandle < handle
    properties
        injectors
        producers
    end
    properties (Dependent, SetObservable)
        Visible              % 'on' / 'off'
        visibleInjectors     % logic array   
        visibleProducers
        Parent
    end
    
    methods
        
        function s = WellPlotHandle(G, W, varargin)
            opt = struct('Parent',     [], ...
                         'Visible',  'on');
            [opt, extraOpt] = merge_options(opt, varargin{:});
              
            if isempty(opt.Parent)
                opt.Parent = gca;
            end
            
            if ~isfield(W(1), 'trajectory')
                W = addTrajectories(W, G, 10);
            end
            %wp = plotWellData(G, W, extraOpt{:}, 'linePlot', true);
            %wp = setParents(wp, opt.Parent);
            %wp = wp(:);
            
            %s.injectors = wp([W.sign]>0);
            %s.producers = wp([W.sign]<0);
            
            s.injectors = plotWellData(G, W([W.sign]>0), extraOpt{:}, 'linePlot', true, 'color', [.2 .2 .7]);
            s.injectors = setParents(s.injectors(:), opt.Parent);
            s.producers = plotWellData(G, W([W.sign]<0), extraOpt{:}, 'linePlot', true, 'color', [.7 .2 .2]);
            s.producers = setParents(s.producers(:), opt.Parent);
            s.Visible   = opt.Visible; 
            set([s.injectors.label], 'interpreter', 'none')
            set([s.producers.label], 'interpreter', 'none')
        end
        
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
        
        function set.visibleInjectors(s, val)
            if ischar(val)
                if strcmp(val, 'on')
                    val = true(numel(s.injectors), 1);
                else
                    val = false(numel(s.injectors), 1);
                end
            elseif ~islogical(val)
                tmp = false(numel(s.injectors), 1);
                tmp(val) = true;
                val = tmp;
            end
            ix = val~=s.visibleInjectors;
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
        end
        function val = get.visibleInjectors(s)
            val = arrayfun(@(x)strcmp(x.label.Visible, 'on'), s.injectors);
            val = val(:);
        end
            
        function set.visibleProducers(s, val)
            if ischar(val)
                if strcmp(val, 'on')
                    val = true(numel(s.producers), 1);
                else
                    val = false(numel(s.producers), 1);
                end
            elseif ~islogical(val)
                tmp = false(numel(s.producers), 1);
                tmp(val) = true;
                val = tmp;
            end
            ix = val~=s.visibleProducers;
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
        end
        function val = get.visibleProducers(s)
            val = arrayfun(@(x)strcmp(x.label.Visible, 'on'), s.producers);
            val = val(:);
        end       
    end
end

function wp = setParents(wp, p)
fn = fieldnames(wp);
for k = 1:numel(fn)
    set([wp.(fn{k})], 'Parent', p)
end
end
    
function W = addTrajectories(W, G, np)
% add well trajectory for plotting, currently just np points along a single
% quadratic curve

for k = 1:numel(W)
    c  = G.cells.centroids(W(k).cells, :);
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

%{
Copyright 2009-2020 SINTEF Digital, Mathematics & Cybernetics.

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
