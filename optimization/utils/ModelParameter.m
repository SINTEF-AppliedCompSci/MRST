classdef ModelParameter
    properties
        name
        type          = 'value';
        boxLims
        distribution  = 'cell';
        Indx
        scaling       = 'linear'
    end
    
    methods
        function p = ModelParameter(model, W, varargin)
            p = merge_options(p, varargin{:});
            p = setupDefaults(p, model, W);
        end
        
        function v = scale(p, pval)
            if strcmp(p.type, 'value') && strcmp(p.scaling, 'linear')
                v = (pval-p.boxLims(1))./diff(p.boxLims, [], 2);
            else
                error('Not implemented')
            end
        end
        
        function pval = unscale(p, v)
            if strcmp(p.type, 'value') && strcmp(p.scaling, 'linear')
                pval = v.*diff(p.boxLims, [], 2) + p.boxLims(1);
            else
                error('Not implemented')
            end
        end
        
        function u = convertToOptimVector(p, pval)
            % Convert parameter pv in model to control vector
            np = numel(pval);
            u = cell(np,1);
            for k = 1:np
                switch p.distribution
                    case 'cell' %parameter distribution per cell
                        switch p.name
                            case {'porevolume','initSw','transmissibility'}
                                u = p.scale(pval);
                            otherwise
                                error('Parameter %s is not implemented',p.name)
                        end
                    case  'connection'
                        switch p.name
                            case {'transmissibility','porevolume','permeability','conntrans'}
                                u = p.scale(pval);
                            otherwise
                                error('Parameter %s is not implemented',p.name)
                        end
                    case 'general'
                        switch p.name
                            case {'transmissibility','porevolume','permeability',...
                                    'swl','swcr', 'swu', 'sgl', ...
                                    'sgcr','sgu','sowcr','sogcr',...
                                    'krw','kro','krg'}
                                u = p.scale(pval);
                            case 'conntrans'
                                u = p.scale(pval);
                            otherwise
                                error('Parameter %s is not implemented',p.name)
                        end
                    otherwise
                        error('Parameter distribution %s is not implemented',p.distribution)
                end
            end
        end
    end
end

function p = setupDefaults(p, model, W)
range = @(x)[.5*min(min(x)), 2*max(max(x))];
switch p.name
    case 'transmissibility'
        if isempty(p.boxLims)
            p.boxLims = range(model.operators.T);
        end
        if isempty(p.Indx)
            p.Indx = (1:numel(model.operators.T))';
        end
    case 'permeability'
        if isempty(p.boxLims)
            p.boxLims = range(model.rock.perm);
        end
        if isempty(p.Indx)
            p.Indx = (1:model.G.cells.num)';
        end
    case 'porevolume'
        if isempty(p.boxLims)
            p.boxLims = range(model.operators.pv);
        end
        if isempty(p.Indx)
            p.Indx = (1:model.G.cells.num)';
        end
    case 'conntrans'
        p.distribution = 'general';
        if isempty(p.boxLims)
            nconn = arrayfun(@(w)numel(w.cells), W);
            tmp = applyFunction(@(x)range(x.WI), W);
            p.boxLims = rldecode(vertcat(tmp{:}), nconn);
        end
        if isempty(p.Indx)
            nconn = arrayfun(@(w)numel(w.cells), W);
            wno   = rldecode((1:numel(W))', nconn(:));
            cno   = mcolon(ones(numel(W),1), nconn);
            p.Indx = [wno(:), cno(:)];
        end
    otherwise
        if isempty(p.boxLims) || isempty(p.Indx)
            error('Can''t set default properties for uknown parameter: %s', p.name);
        end
end
end


