classdef ModelParameter
    properties
        name
        type          = 'value';    % 'value'/'multiplier'
        boxLims                     % upper/lower value(s) for parameters (used for scaling)
        subset                      % subset of parameters (or subset of wells)
        scaling       = 'linear'    % 'linear'/'log'
        referenceValue              % parameter reference values (used for type 'multiplier') 
        belongsTo                   % model/shcedule/state0/
        location                    % e.g., {'model', 'operators', 'T'}
        nParam                      % number of parameters
        lumping                     % parameter lumping vector (partition vector) 
        setfun                      % possibly custom set-function (default is setfield)
        scalingBase = nan;          
        controlSteps = [];          % Default set/get for all control steps
    end
    
    methods
        function p = ModelParameter(setup, varargin)
            [p, extra] = merge_options(p, varargin{:});
            assert(~isempty(p.name), 'Parameter name can''t be defaulted');
            if isempty(p.belongsTo) || isempty(p.location) 
                p = setupByName(p, setup);
            end
            opt = struct('relativeLimits', [.5 2]);
            opt = merge_options(opt, extra{:});
            p   = setupDefaults(p, setup, opt);
            if isempty(p.setfun)
                % use default
                p.setfun = @(obj, loc, v)setfield(obj, loc{:}, v);
            end
            checkSetup(p, setup);
        end
        %------------------------------------------------------------------
        function vs = scale(p, pval)
            % map parameter pval to "control"-vector v \in [0,1]
            vs = (pval-p.boxLims(:,1))./diff(p.boxLims, [], 2);
            if strcmp(p.scaling, 'exp')
                vs = expScale(vs, p);
            elseif strcmp(p.scaling, 'log')
                vs = logScale(vs, p);
            end
        end
        %------------------------------------------------------------------
        function pval = unscale(p, vs)
            % retrieve parameter pval from "control"-vector v \in [0,1]
            pval = vs;
            if strcmp(p.scaling, 'exp')
                pval = logScale(pval, p);
            elseif strcmp(p.scaling, 'log')
                pval = expScale(pval, p);
            end
            pval = pval.*diff(p.boxLims, [], 2) + p.boxLims(:,1);
        end
        %------------------------------------------------------------------
        function gs = scaleGradient(p, g, pval)
            % map gradient wrt param to gradient vs "control"-vector
            % parameter value pval only needed for log-scalings
            gs = g.*diff(p.boxLims, [], 2);
            if ~strcmp(p.scaling, 'linear')
                tmp = (pval-p.boxLims(:,1))./diff(p.boxLims, [], 2);
                if strcmp(p.scaling, 'exp')
                    gs = gs./dExpScale(tmp, p);
                elseif(strcmp(p.scaling, 'log'))
                    gs = gs./dLogScale(tmp, p);
                end
            end
        end
        %------------------------------------------------------------------
        function g = collapseGradient(p, g)
            assert(strcmp(p.belongsTo, 'state0'),['This function is only',...
                'intended to collapse gradient arrays from parameters that',...
                'belongs to state0'])
            % take sum of each lump
            if ~isempty(p.lumping) && isnumeric(p.lumping)
                g = accumarray(p.lumping, g(p.subset), [], @sum); 
            end
        end
        %------------------------------------------------------------------
        function v = getParameter(p, setup)
            if ~strcmp(p.type, 'multiplier')
                v = p.getParameterValue(setup);
            else
                v = p.getParameterValue(setup, false)./p.referenceValue;
                v = collapseLumps(v, p.lumping);
            end
        end
        %------------------------------------------------------------------
        function setup = setParameter(p, setup, v)
            if ~strcmp(p.type, 'multiplier')
                setup = p.setParameterValue(setup, v);
            else
                v  = expandLumps(v, p.lumping).*p.referenceValue;
                setup = p.setParameterValue(setup, v, false);
            end
        end
        
        %------------------------------------------------------------------
        function v = getParameterValue(p, setup, doCollapse)
            if nargin < 3
                doCollapse = true;
            end
            if ~(strcmp(p.belongsTo, 'schedule') && strcmp(p.location{1}, 'control'))
                v = getfield(setup.(p.belongsTo), p.location{:});
                v = v(p.subset);
                if doCollapse
                    v = collapseLumps(v, p.lumping);
                end
            else % control-parameter (assume constant selection of control steps)
                control = setup.schedule.control(p.controlSteps(1));
                v  = p.getControlParameterValue(control, doCollapse);
            end
        end
        %------------------------------------------------------------------       
        function setup = setParameterValue(p, setup, v, doExpand)
            if nargin < 4
                doExpand = true;
            end
            if doExpand
                v  = expandLumps(v, p.lumping);
            end
            if ~(strcmp(p.belongsTo, 'schedule') && strcmp(p.location{1}, 'control'))
                if isnumeric(p.subset)
                    tmp = getfield(setup.(p.belongsTo), p.location{:});
                    v   = setSubset(tmp, v, p.subset);
                end
                setup.(p.belongsTo) = ...
                    p.setfun(setup.(p.belongsTo), p.location, v);
            else % control-parameter (assume constant over selected control steps)
                nc = numel(p.controlSteps);
                for k = 1:nc
                    step = p.controlSteps(k);
                    setup.schedule.control(step) = ...
                        p.setControlParameterValue(setup.schedule.control(step), v);
                end
            end
        end
        %------------------------------------------------------------------       
        function v = getControlParameterValue(p, control, doCollapse)
            if nargin < 3
                doCollapse = true;
            end
            if strcmp(p.location{2}, 'W') 
                % well-parameter subset applies to well numbers 
                loc = p.location(3:end);
                v = applyFunction(@(x)getfield(x, loc{:}), control.W(p.subset));
                v = vertcat(v{:});
            else
                % apply subset to parameter
                loc = p.location(2:end);
                v = getfield(control, loc{:});
                v = v(p.subset);
            end      
            if doCollapse
                v = collapseLumps(v, p.lumping);
            end
        end
        %------------------------------------------------------------------
        function control = setControlParameterValue(p, control, v)
            if strcmp(p.location{2}, 'W')
                % well-parameter special treatment
                sub = p.subset;
                if ~isnumeric(sub)
                    sub = 1:numel(control.W);
                end
                nc = arrayfun(@(w)numel(w.cells), control.W(sub));
                [i1, i2] = deal(cumsum([1;nc(1:end-1)]), cumsum(nc));
                v  = applyFunction(@(i1,i2)v(i1:i2), i1, i2);
                loc = p.location(3:end);
                for k = 1:numel(sub)
                    control.W(sub(k)) = p.setfun(control.W(sub(k)), loc, v{k});
                end
            else
                % other parameter (set subset of parameter)
                loc = p.location(2:end);
                if isnumeric(p.subset)
                    tmp = getfield(control, loc{:});
                    v   = setSubset(tmp, v, p.subset);
                end
                control = p.setfun(control, loc, v);
            end
        end
    end
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function p = setupDefaults(p, setup, opt)
% Make sure setup makes sense and add boxLims if not provided
if isempty(p.subset)
    p.subset = ':';
end
if islogical(p.subset)
    p.subset = find(p.subset);
end

if strcmp(p.belongsTo, 'schedule')
    if isempty(p.controlSteps)
        p.controlSteps = (1:numel(setup.schedule.control));
    else
        % set name depending on first selected control step
        p.name = sprintf('%s_step%d', p.name, p.controlSteps(1));
    end
end

v = getParameterValue(p, setup, false);

if ~isempty(p.lumping)
    p.nParam = max(p.lumping);
else
    p.nParam = numel(v);
end
   
if strcmp(p.type, 'multiplier')
    p.referenceValue = v;
end

range = @(x)[min(min(x)), max(max(x))];
if isempty(p.boxLims)
    rlim  = opt.relativeLimits;
    if strcmp(p.type, 'value')
        p.boxLims = range(v).*rlim;
    else
        p.boxLims = rlim;
    end
    % special treatment of saturations
    if any(strcmp(p.name, {'sw', 'sg'}))
         p.boxLims = [0, 1];
    end
end

assert(any(size(p.boxLims,1) == [1, p.nParam]), ...
    'Property ''boxLims'' does not match number of parameters');
end

%--------------------------------------------------------------------------
function p = setupByName(p, SimulatorSetup)
% setup for typical parameters
setfun = [];
switch lower(p.name)
    case 'transmissibility'
        belongsTo = 'model';
        location = {'operators', 'T'};
    case {'permx', 'permy', 'permz'}
        belongsTo = 'model';
        col = find(strcmpi(p.name(end), {'x', 'y', 'z'}));
        location = {'rock', 'perm', {':', col}};
        setfun   = @setPermeabilityFun;
    case 'porevolume'
        belongsTo = 'model';
        location = {'operators', 'pv'};
    case 'conntrans'
        belongsTo = 'schedule';
        location = {'control', 'W', 'WI'};
    case {'sw', 'sg'}
        belongsTo = 'state0';
        col = SimulatorSetup.model.getPhaseIndex(upper(p.name(end)));
        location = {'s', {':', col}};
        oix = SimulatorSetup.model.getPhaseIndex('O');
        assert(~isempty(oix), ...
            'Current assumption is that oil is the dependent phase');
        setfun   = @(obj, loc, v)setSaturationFun(obj, loc, v, oix);
    case 'pressure'
        belongsTo = 'state0';
        location = {'pressure'};
    case {'swl', 'swcr', 'swu', 'sowcr', 'sogcr', 'sgl', 'sgcr', ...
            'sgu', 'krw', 'kro', 'krg'}
        belongsTo = 'model';
        map = getScalerMap();
        ix  = map.kw.(upper(p.name));
        [ph, col] = deal(map.ph{ix(1)}, ix(2));
        location = {'rock', 'krscale', 'drainage', ph, {':', col}};
        setfun   = @setRelPermScalersFun;
    otherwise
        error('No default setup for parameter: %s\n', p.name);     
end
if isempty(p.belongsTo)
    p.belongsTo = belongsTo;
end
if isempty(p.location)
    p.location = location;
end
if isempty(p.setfun) && ~isempty(setfun)
    p.setfun = setfun;
end
end

%--------------------------------------------------------------------------            
function map = getScalerMap()
phOpts = {'w', 'ow', 'g', 'og'};
kw  = struct('SWL',   [1,1], 'SWCR',  [1,2], 'SWU', [1,3], ...
             'SGL',   [3,1], 'SGCR',  [3,2], 'SGU', [3,3], ...
             'SOWCR', [2,2], 'SOGCR', [4,2], ...
             'KRW',   [1,4], 'KRO',   [2,4], 'KRG', [3,4]);
map = struct('ph', {phOpts}, 'kw', kw);
end

%--------------------------------------------------------------------------
function v = collapseLumps(v, lumps)
% take mean of each lump
if ~isempty(lumps) && isnumeric(lumps)
    if isa(v, 'double')
        v = accumarray(lumps, v, [], @mean);
    else % special treatment in case of ADI
        M = sparse(lumps, (1:numel(lumps))', 1);
        v = (M*v)./sum(M,2);
    end
end
end

%--------------------------------------------------------------------------
function v = expandLumps(v, lumps)
if ~isempty(lumps) && isnumeric(lumps)
    v = v(lumps);
end
end

%--------------------------------------------------------------------------
function v = setSubset(v, vi, sub)
if isa(vi, 'ADI')
    v = double2ADI(v, vi);
end
v(sub) = vi;
end

%--------------------------------------------------------------------------
function checkSetup(p, setup)
if ~isempty(p.lumping)
     tmp = p.getParameterValue(setup, false);
     assert(numel(p.lumping) == numel(tmp), 'Lumping vector has incorrect size');
end
tmp = p.getParameter(setup);
assert(numel(tmp)==p.nParam, 'Report error to develolper')
assert(any(size(p.boxLims,1) == [1, p.nParam]), ...
       'The number of upper/lower limits should be 1 or nParam');
assert(all(tmp >= p.boxLims(:,1) & tmp <= p.boxLims(:,2)), ...
       'Parameter values are not within given limits')
end

%--------------------------------------------------------------------------    
function model = setPermeabilityFun(model, location, v)
% utility for setting permx/y/z possibly as AD and include effect on
% transmissibilities
[nc, nd] = size(model.rock.perm);
perm = model.rock.perm;
if ~iscell(perm)
    perm = mat2cell(perm, nc, ones(1,nd));
end
col = location{end}{end};
assert(col<=nd, 'Can''t get column %d since perm has %d column(s)', col, nd);
perm{col} = v;
% transmissibilities
th = 0;
for k = 1:nd
    th = th + perm2directionalTrans(model, perm{k}, k);
end
cf = model.G.cells.faces(:,1);
nf = model.G.faces.num;
% mapping from from cell-face to face
M = sparse(cf, (1:numel(cf))', 1, nf, numel(cf));
% consider only internal faces
ie = model.operators.internalConn;
model.operators.T = 1./(M(ie,:)*(1./th));
% if all perms are doubles, set to nc x nd array
if all(cellfun(@(p)isa(p, 'double'), perm))
    model.rock.perm = cell2mat(perm);
else
    model.rock.perm = perm;
end
end

%--------------------------------------------------------------------------    
function state = setSaturationFun(state, location, v, oix)
assert(isa(v, 'double'), 'Setting saturation to class %s is not supported', class(v));
pix = location{end}{end};
ds = v-state.s(:, pix);
state.s(:, pix) = v;
state.s(:, oix) =  state.s(:, oix) - ds;
end

%--------------------------------------------------------------------------   
function model = setRelPermScalersFun(model, location, v)
if ~isa(v, 'ADI')
    model = setfield(model, location{:}, v);
else
    % last location is column no, second last is phase
    col = location{end}{end};
    ph  = location{end-1};
    d   = getfield(model, location{1:end-2});  %#ok
    if ~isfield(d, 'tmp') || ~isfield(d.tmp, ph)
        nc = model.G.cells.num;
        d.tmp.(ph) = mat2cell(d.(ph), nc, ones(1,4));
    end
    d.tmp.(ph){col} = v;
    d.(ph) = @(cells, col)d.tmp.(ph){col}(cells);
    model.rock.krscale.drainage = d;
end
end

%--------------------------------------------------------------------------          
function ti = perm2directionalTrans(model, p, cdir)
% special utility function for calculating transmissibility along coordinate direction cdir
% In particular:
% trans = t1+t2+t2, where ti = perm2dirtrans(model, perm(:, i), i);
assert(size(p,2)==1, ...
       'Input p should be single column representing permeability in direction cdir');
% make perm represent diag perm tensor with zero perm orthogonal to cdir
dp = value(p);
r.perm = zeros(numel(dp), 3);
r.perm(:, cdir) = dp;
ti = computeTrans(model.G, r);
if isa(p, 'ADI')
    % make ti ADI (note ti is linear function of p)
    p = p./dp;
    cellno = rldecode(1 : model.G.cells.num, diff(model.G.cells.facePos), 2).';
    ti = ti.*p(cellno);
end
end

%--------------------------------------------------------------------------       
function vs = expScale(v, p)
base = getScalingBase(p);
vs = (base.^v - 1)./(base - 1);
end

%--------------------------------------------------------------------------       
function ds = dExpScale(v, p)
base = getScalingBase(p);
ds = (base.^v).*(log(base)./(base - 1));
end

%--------------------------------------------------------------------------       
function vs = logScale(v, p)
base = getScalingBase(p);
vs = log((base-1).*v+1)./log(base);
end

%--------------------------------------------------------------------------       
function ds = dLogScale(v, p)
base = getScalingBase(p);
ds = (base-1)./( ((base-1).*v+1).*log(base) );
end

%--------------------------------------------------------------------------     
function base = getScalingBase(p)
base = p.scalingBase;
if isnan(base)
    base = p.boxLims(:,2)./p.boxLims(:,1);
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
