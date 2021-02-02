classdef ModelParameter
    properties
        name
        type          = 'value'; % 'value'/'multiplier'
        boxLims                  % upper/lower value(s) for parameters
        distribution  = 'cell';  % not sure we need this one
        Indx                     % 
        scaling       = 'linear' % 'linear'/'log'
        initialValue             % only used for multipliers
        belongsTo                % model/well/state0
        location                 % e.g., {'operators', ''}
        n                        % number of 
        lumping                  % 
    end
    
    methods
        function p = ModelParameter(problem, varargin)
            [p, extra] = merge_options(p, varargin{:});
            assert(~isempty(p.name), 'Parameter name can''t be defaulted');
            [p.belongsTo, p.location] = setupAddress(p.name);
            opt = struct('relativeLimits', [.5 2]);
            opt = merge_options(opt, extra{:});
            p = setupDefaults(p, problem, opt);
        end
        
        function vs = scale(p, pval)
            % map parameter pval to "control"-vector v \in [0,1]
            if strcmp(p.type, 'multiplier')
                pval = pval./p.initialValue;
            end
            pval = collapseLumps(p, pval, @mean);
            if strcmp(p.scaling, 'linear')
                vs = (pval-p.boxLims(:,1))./diff(p.boxLims, [], 2);
            elseif strcmp(p.scaling, 'log')
                logLims = log(p.boxLims);
                vs = (log(pval)-logLims(:,1))./diff(logLims, [], 2);
            end
        end
     
        function pval = unscale(p, vs)
            % retrieve parameter pval from "control"-vector v \in [0,1]
            if strcmp(p.scaling, 'linear')
                pval = vs.*diff(p.boxLims, [], 2) + p.boxLims(:,1);
            elseif strcmp(p.scaling, 'log')
                logLims = log(p.boxLims);
                pval = exp(vs.*diff(logLims, [], 2) + logLims(:,1));
            end
            pval = expandLumps(p,pval);
            if strcmp(p.type, 'multiplier')
                pval = pval.*p.initialValue;
            end
        end
        
        function gs = scaleGradient(p, g, pval)
            % map gradient wrt param to gradient vs "control"-vector
            % parameter value pval only needed for log-scalings
            g = collapseLumps(p, g, @sum);
            if strcmp(p.scaling, 'linear')
                if strcmp(p.type, 'value')
                    gs = g.*diff(p.boxLims, [], 2);
                elseif strcmp(p.type, 'multiplier')
                    gs = (g.*p.initialValue).*diff(p.boxLims, [], 2);
                end
            elseif strcmp(p.scaling, 'log')
                gs = (g.*pval).*log(diff(p.boxLims, [], 2));
            end
        end
        
        function v = getParameterValue(p, problem)
            if strcmp(p.belongsTo, 'model')
                v = p.getModelParameterValue(problem.model);
            elseif strcmp(p.belongsTo, 'well')
                v = p.getWellParameterValue(problem.schedule.control(1).W);
            end
        end
                
        function problem = setParameterValue(p, problem, v)
            if strcmp(p.belongsTo, 'model')
                problem.model = p.setModelParameterValue(problem.model, v);
            elseif strcmp(p.belongsTo, 'well')
                for k = 1:numel(problem.schedule.control)
                    problem.schedule.control(k).W = ...
                        p.setWellParameterValue(problem.schedule.control(k).W, v);
                end
            end
        end
        
        function v = getModelParameterValue(p, model)
            assert(strcmp(p.belongsTo, 'model'))
            v = getfield(model, p.location{:});
            if ~isempty(p.Indx) || ~isnumeric(p.Indx)
                v = v(p.Indx);
            end
        end
        
        function model = setModelParameterValue(p, model, v)
            assert(strcmp(p.belongsTo, 'model'))
            if ~isempty(p.Indx) || ~isnumeric(p.Indx)
                model = setfield(model, p.location{:}, v);
            else
                tmp = getfield(model, p.location{:});
                tmp(p.Indx) = v;
                model = setfield(model, p.location{:}, tmp);
            end
        end
        
        function v = getWellParameterValue(p, W, cellOutput)
            assert(strcmp(p.belongsTo, 'well'))
            v = applyFunction(@(x)getfield(x, p.location{:}), W);
            if ~(nargout == 3) || ~cellOutput
                v = vertcat(v{:});
            end
        end
        
        function W = setWellParameterValue(p, W, v)
            if ~iscell(v)
                nc = arrayfun(@(w)numel(w.cells), W);
                v = mat2cell(v, nc, 1);
            end
            for k = 1:numel(W)
                W(k) = setfield(W(k), p.location{:}, v{k});
            end
        end
        
        function m = getMultiplerValue(p, problem, doLump)
            if strcmp(p.type, 'multiplier')
                m = p.getParameterValue(problem)./p.initialValue;
                if nargin == 3 && doLump
                    m = collapseLumps(p, m, @mean);
                end
            else
                error('Parameter %s is not of type ''multiplier''', p.name);
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

function p = setupDefaults(p, problem, opt)
rlim = opt.relativeLimits;
range = @(x)[min(min(x)), max(max(x))];
if strcmp(p.belongsTo, 'model')
    v = p.getModelParameterValue(problem.model);
    if isempty(p.boxLims)
        if strcmp(p.type, 'value')
            p.boxLims = range(v).*rlim;
        else
            p.boxLims = rlim;
        end
    end
    if isempty(p.Indx)
        p.Indx = ':';
    end
    if strcmp(p.type, 'multiplier')
        p.initialValue = v;
    end
    p.n = numel(v);
elseif strcmp(p.belongsTo, 'well')
    % we only have WI at the moment so this is a special case
    assert(strcmp(p.name, 'conntrans'), 'Needs updating for: %s', p.name);
    p.distribution = 'general';
    W = problem.schedule.control(1).W;
    nconn = arrayfun(@(w)numel(w.cells), W);
    if isempty(p.boxLims)
        nconn = arrayfun(@(w)numel(w.cells), W);
        if strcmp(p.type, 'value')
            tmp = applyFunction(@(x)range(x.WI).*rlim, W);
            p.boxLims = rldecode(vertcat(tmp{:}), nconn);
        else
            p.boxLims = repmat(rlim, sum(nconn), 1);
        end
    end
    if isempty(p.Indx)
        wno   = rldecode((1:numel(W))', nconn(:));
        cno   = mcolon(ones(numel(W),1), nconn);
        p.Indx = [wno(:), cno(:)];
    end
    if strcmp(p.type, 'multiplier')
        p.initalValue = vertcat(W.WI);
    end
    p.n = sum(nconn);
end
end


function [belongsTo, location] = setupAddress(name)
    belongsTo = 'model';
    switch lower(name)
        case 'transmissibility'
            location = {'operators', 'T'};
        case 'permeability'
            location = {'rock', 'perm'};
        case 'porevolume'
            location = {'operators', 'pv'};
        case 'conntrans'
            belongsTo = 'well';
            location = {'WI'};
        case {'swl', 'swcr', 'swu', 'sowcr', 'sogcr', 'sgl', 'sgcr', ...
              'sgu', 'krw', 'kro', 'krg'}
            map = getScalerMap();
            ix  = map.kw.(upper(name));
            [ph, col] = deal(map.ph{ix(1)}, ix(2));
            location = {'rock', 'krscale', 'drainage', ph, {':', col}};
        otherwise
            error('Uknown parameter name: %s', name);
    end
end
            
function map = getScalerMap()
phOpts = {'w', 'ow', 'g', 'og'};
kw  = struct('SWL',   [1,1], 'SWCR',  [1,2], 'SWU', [1,3], ...
             'SGL',   [3,1], 'SGCR',  [3,2], 'SGU', [3,3], ...
             'SOWCR', [2,2], 'SOGCR', [4,2], ...
             'KRW',   [1,4], 'KRO',   [2,4], 'KRG', [3,4]);
map = struct('ph', {phOpts}, 'kw', kw);
end

function v = collapseLumps(p, v, op)
% do ommitnan for safety
fn = @(x)op(x, 'omitnan');
if ~isempty(p.lumping) || ~isnumeric(p.lumping)
    v = accumarray(p.lumping, [], fn);
end
end

function v = expandLumps(p, v)
if ~isempty(p.lumping) || ~isnumeric(p.lumping)
    v = v(p.lumping);
end
end
    
