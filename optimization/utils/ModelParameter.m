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
        grad_location            % e.g., {'init',''}
        n                        % number of 
        lumping                  %
        initial_val
    end
    
    methods
        function p = ModelParameter(problem, varargin)
            [p, extra] = merge_options(p, varargin{:});
            assert(~isempty(p.name), 'Parameter name can''t be defaulted');
            [p.belongsTo, p.location, p.grad_location] = setupAddress(p.name);
            opt = struct('relativeLimits', [.5 2]);
            opt = merge_options(opt, extra{:});
            p = setupDefaults(p, problem, opt);
        end
        
        function vs = scale(p, pval)
            % map parameter pval to "control"-vector v \in [0,1]
            if strcmp(p.type, 'multiplier')
                pval = pval./p.initialValue;
            end
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
            elseif strcmp(p.belongsTo, 'state0')
                v = p.getState0ParameterValue(problem.state0);
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
            elseif strcmp(p.belongsTo, 'state0')
                problem.state0 = p.setState0ParameterValue(problem.state0, v);
            end
        end
        
        function v = getModelParameterValue(p, model)
            assert(strcmp(p.belongsTo, 'model'))
            v = getfield(model, p.location{:});
            v = collapseLumps(p, v, @mean);
        end
        
        function model = setModelParameterValue(p, model, v)
            assert(strcmp(p.belongsTo, 'model'))
            tmp = getfield(model, p.location{:});            
            v = expandLumps(p,v,tmp); 
                model = setfield(model, p.location{:}, v);

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
        
        function v = getState0ParameterValue(p, state0)
            assert(strcmp(p.belongsTo, 'state0'))
            v = getfield(state0, p.location{:});            
            v = collapseLumps(p, v(:,1), @mean);
        end
        
        function state0 = setState0ParameterValue(p, state0, v)
            assert(strcmp(p.belongsTo, 'state0'))
            tmp = getfield(state0, p.location{:});            
            v = expandLumps(p,v,tmp(:,1));
            if strcmp(lower(p.name), 'initsw')
                    v = [v,1-v]; %TODO: extend it to more phases
            end
            state0 = setfield(state0, p.location{:},v);
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
                                u = p.scale(pval);                          
            end
        end
    end
end

function p = setupDefaults(p, problem, opt)
rlim = opt.relativeLimits;
range = @(x)[min(min(x)), max(max(x))];
if (strcmp(p.belongsTo, 'model') || strcmp(p.belongsTo, 'state0'))
    v = getParameterValue(p, problem);
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
    if ~isempty(p.lumping) || ~isnumeric(p.lumping)
        p.n = numel(unique(p.lumping(~isnan(p.lumping)))); %Number of unique parameters
    else
        p.n =  numel(v);
    end
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


function [belongsTo, location, grad_location] = setupAddress(name)
    belongsTo = 'model';
    grad_location = {name};
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
        case 'initsw'
            belongsTo = 'state0';
            location = {'s'};
            grad_location = {'init','sW'};
        case 'p0'
            belongsTo = 'state0';
            location = {'pressure'};
            grad_location = {'init','pressure'};
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
% do ommitnan for safety on v
fn = @(x)op(x, 'omitnan');

if ~isempty(p.lumping) || ~isnumeric(p.lumping)
    v = accumarray_nan(p.lumping,v,fn); % lumping may have Nan values
end
end

function model_val = expandLumps(p, v, model_val)
    if ~isempty(p.lumping) || ~isnumeric(p.lumping)        
        for i =  1:length(model_val)
            if ~isnan(p.lumping(i))
                model_val(i,1) = v(p.lumping(i));
            end
        end
    else 
        model_val = v;
    end
end
    
function v = accumarray_nan(subs,val, fn)
    if any(isnan(subs))
       subs(isnan(subs)) = 0;
       for i = 1:length(unique(subs))-1
           Indx = find(subs==i);
           v(i,1) = fn(val(Indx));
       end
    else
       v = accumarray(subs,val,[],fn);
    end
end

