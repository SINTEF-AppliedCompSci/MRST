function varargout = computeEnsembleDiagnostics(setupFn, varargin)
% General purpose routine for computing various diagnostics on ensemble of 
% models.   
%
% SYNOPSIS:
% [states, diagnostics] = computeEnsembleDiagnostics(setupFn, 'pn1', pv1)
% computeEnsembleDiagnostics(setupFn, 'outputDir', '...', 'pn1', pv1)
%
% REQUIRED PARAMETERS:
% setupFn - function handle to setup-function. Required to be of format
%           [model_i, W_i] = setupFn(i), or
%           [model_i, W_i, state0_i] = setupFn(i),
%           where i is ralization number 
%
% OPTIONAL PARAMETERS:
% 'control'   - control struct with fields 'W', 'bc' and/or 'src'. If
%               non-empty diagnostics will be computed with
%               model-independent fields from control and model-dependent
%               fields from W_i (e.g., 'WI'). (Default [])
% 
% 'outputDir' - Output directory. If empty, output is not written to file.
%               If string, output is written to outputDir/realization1, 
%               outputDir/realization2, etc,
%               If function handle, output is written to outputDir(1),
%               outputDir(2), etc. (Default '')
%
% 'fluid'     - Fluid (AD-type). If defaulted, model.fluid will be used
%
% 'state0'    - Initial state. Will only be used if nargout(setupFn) == 2 
%
% 'ellipticSolver'- Linear solver passed to incompTPFA
%
% 'ensembleSize' - Number of realizations. If undefined, setupFn is assumed
%                  to provide empty output for n > ensembleSize. 
%
% 'memberIx'  - indices to ensemble subset. (default 1:ensembleSize)
% 
% 'overWrite' - Overwrite existing data in outputDir without asking. 
%               (Defailt false)
%
% 'rtd'       - options are 'estimate', 'compute', 'both', 'none' (default) 
%
% 'rtdOpts'   - additional options passed to computeRTD/estimateRTD
%
% 'diagnFn'   - function handle or list of handles to additional routines 
%               assumed to be of format 
%                   output = diagnFn(state, model, D, WP).
%               Output as field in diagnostics or/and written as seperate
%               file. Name of field/file is func2str(diagnFn)
%
% 'p_ref'     - Reference pressure for convertion to incompressible fluid.
%
% [additional]- passed to computeTOFandTracer 
% 
%
% RETURNS
% If nargout == 2:
% state         - contining computed pressure/flux
% diagnostics   - structure containing fields 'D' (tof/tracer fields), 'WP'
%                 (interaction allocation/volumes) and 'wellCommuniaction' 
%                 Additional fields according to diagnFn.

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

opt = struct('control',                  [], ...
             'outputDir',                '', ...
             'state0',                   [], ...
             'ellipticSolver',    @mldivide, ...
             'ensembleSize',            inf, ...
             'memberIx',                 [], ...
             'overWrite',             false, ...
             'rtd',              'estimate', ...  %'none', 'estimate', 'compute', 'both' 
             'rtdOpts',                {{}}, ...
             'diagnFn',                  '', ...         
             'maxTOF',                   [], ...
             'computeWellTOFs',        true, ...
             'processCycles',          true, ...
             'firstArrival',           true, ...
             'writeToDiskMode',       'all', ...   % 'none', 'light', 'all'
             'p_ref',             400*barsa, ...
             'showWaitbar',           true);     

[opt, extra_opt] = merge_options(opt, varargin{:}); 
opt.rtdOpts = [opt.rtdOpts, {'showWaitbar', opt.showWaitbar}];

if nargout > 0
    stct = struct('D', [], 'WP', [], 'wellCommunication', [], 'rtd', []);
    if opt.ensembleSize < inf
        states      = cell(opt.ensembleSize, 1);
        diagnostics = repmat(stct, [opt.ensembleSize, 1]);
    else
        states      = cell(1,1);
        diagnostics = stct;
    end
end

% check outputDir
if ~strcmp(opt.writeToDiskMode, 'none')
    outputDir = opt.outputDir;
%     if ~opt.overWrite
%         if isfolder(outputDir) && ~isempty(dir(fullfile(outputDir, 'D*.mat')))
%             answ = questdlg('Output folder is non-empty', 'Warning', ...
%                 'Proceed anyway', 'Cancel', 'Cancel');
%             if strcmp(answ, 'Cancel')
%                  return;
%              end
%         end
%     end
    if strcmp(opt.writeToDiskMode, 'light')
       h = getHandlers(outputDir, {'WP', 'RTD', 'other','statistics'});
    else
       h = getHandlers(outputDir, {'states', 'D', 'WP', 'RTD', 'other','statistics'}); 
    end
end



% handle additional output functions
funcs = struct('handles', {{}}, 'names', {{}});
if ~isempty(opt.diagnFn)
    if ~iscell(opt.diagnFn), opt.diagnFn = {opt.diagnFn}; end
    funcs.handles = opt.diagnFn;
    funcs.names   = cellfun(@func2str, opt.diagnFn, 'UniformOutput', false);
end

memberIx = opt.memberIx;
if isempty(memberIx)
    nMem = min(opt.ensembleSize, 1000); % just something big
    memberIx = (1:nMem)';
end

% check for previously computed diagnostics
tmp = zeros(max(memberIx), 1);
flds = fieldnames(h);
for k = 1:numel(flds)
    ix = h.(flds{k}).getValidIds;
    if max(ix) > numel(tmp)
        tmp(max(ix)) = 0;
    end
    tmp(ix) = tmp(ix)+1;
end
isComputed = (tmp == numel(flds));

% main loop
if opt.showWaitbar
    hwb  = waitbar(0,'Computing individual ensemble members...');
    hwbs = waitbar(0,'Progress for each member');
else
    hwbs = [];
end

        statistics.poro.globalmin = 0;
        statistics.poro.globalmax = 1;

for k = 1:numel(memberIx)
    ix = memberIx(k);
    if isComputed(ix) 
        dispif(mrstVerbose, 'Diagnostics for member no: %d has allready been computed, continuing ...\n', ix);
        continue;
    end
    % setup
    tmp = setupFn(ix);
    if isempty(tmp.model), break; end
    if ~isfield(tmp, 'state0') || isempty(tmp.state0)
        dispif(mrstVerbose, 'Initial state not provided, using default');
        state0 = initResSol(tmp.model.G, opt.p_ref);
    else
        state0 = tmp.state0;
    end
    [model, W] = deal(tmp.model, tmp.W);
    
    % use intput state if not empty
    if ~isempty(opt.state0)
        state0 = opt.state0;
    end
    
    % convert to incompressible fluid
    model.fluid = convertToIncompFluid(model, 'state', state0);
    
    % handle well controls
    [W, bc, src] = getControl(opt.control, W);   
    
    % convert to reservoir-rates for non-unit b-factors
    W = adjustControlsForNonUnitBFactors(W, model, opt.p_ref);
           
    % compute pressure
    state = incompTPFA(state0, model.G, model.operators.T_all, model.fluid, ...
                       'wells', W, 'bc', bc, 'src', src, ...
                       'LinSolve', opt.ellipticSolver, 'use_trans', true);
   
    % make sure flux is column vector
    state.flux = state.flux(:);
                   
    % diagnostics part               
    D  = computeTOFandTracer(state, model.G, model.rock, 'wells', W, ...
                             'maxTOF',            opt.maxTOF, ...
                             'computeWellTOFs',   opt.computeWellTOFs, ...
                             'processCycles',     opt.processCycles, ...
                             'firstArrival',      opt.firstArrival, ...
                             extra_opt{:});
                   
    WP = computeWellPairs(state, model.G, model.rock, W, D);
    salloc = cellfun(@sum, {WP.inj.alloc}, 'UniformOutput',false);
    wellCommunication = vertcat(salloc{:});

    
    % rtd options
    rtd = [];
    if strcmp(opt.rtd, 'estimate')
        rtd = estimateRTD(model.operators.pv, D, WP, opt.rtdOpts{:});
    elseif strcmp(opt.rtd, 'compute') 
        rtd = computeRTD(state, model.G, model.operators.pv, D, WP, W, opt.rtdOpts{:}, 'wbHandle', hwbs);
    elseif strcmp(opt.rtd, 'both')
        if ~isempty(opt.rtdOpts)
            warning('Can''t apply rtd-options both to computeRTD and estimateRTD');
        end
        rtd = {estimateRTD(model.operators.pv, D, WP), computeRTD(state, model.G, model.operators.pv, D, WP)};
    end
        
    

    
    % histogram calculation static rock data
     if (k==1) 
     % Porosity       
        hist_rock.poro = computeHistogram(model.rock.poro,...
                                          'N_bins', 500,...
                                          'BinLimits',[0 1]',...
                                          'Normalization', 'cumcount');
     % Permeability 
        hist_rock.perm = computeHistogram(model.rock.perm,...
                                          'N_bins', 500,...
                                          'Normalization', 'cumcount');
        

        
        statistics.perm.globalmin = hist_rock.perm.min;
        statistics.perm.globalmax = hist_rock.perm.max;
     else  
                                      
        hist_rock.poro = computeHistogram(model.rock.poro,...
                                          'N_bins', 500,...
                                          'BinLimits',[0 1]',...
                                          'Normalization', 'cumcount'); 
        hist_rock.perm = computeHistogram(model.rock.perm,...
                                          'N_bins', 500,...
                                          'Normalization', 'cumcount');


        
        for i=1:size(hist_rock.perm.min,2)
          if hist_rock.perm.min(1,i) > 0
                  statistics.perm.globalmin(1,i) = min(hist_rock.perm.min(1,i), statistics.perm.globalmin(1,i));
          end
          statistics.perm.globalmax(1,i) = max(hist_rock.perm.max(1,i), statistics.perm.globalmax(1,i));
        end
            
     end
     
    % additional diagnostic functions
    output = cell(numel(funcs.handles));
    for j = 1:numel(output)
        output{j} = funcs.handles{j}(state, model, D, WP,statistics);
    end
    
    if nargout > 0
        states{k}      = state;
        diagnostics(k).D  = D;
        diagnostics(k).WP = WP;
        diagnostics(k).wellCommunication = wellCommunication;
        if ~isempty(rtd)
            diagnostics(k).rtd = rtd;
        end
        for j = 1:numel(output)
            diagnostics(k).(funcs.names{j}) = output{j};
        end
    end
     

    if ~strcmp(opt.writeToDiskMode, 'none')
        if strcmp(opt.writeToDiskMode, 'all')
            h.states{ix} = {state};
            h.D{ix}      = {D};
        end
        h.WP{ix} = {WP};
        h.RTD{ix} = {rtd};
        other = struct('wellCommunication', wellCommunication, 'fluid', model.fluid, 'rock',hist_rock);
        for j = 1:numel(output)
            other.(funcs.names{j}) = output{j};
        end
        h.other{ix} = {other};
    end
    if opt.showWaitbar
        waitbar(k/numel(memberIx),hwb);
    end
end
        h.statistics{1} = {statistics};

if opt.showWaitbar
    close([hwb hwbs]);
end
if nargout > 0
    varargout = {states, diagnostics};
end
end

%--------------------------------------------------------------------------

function [W_out, bc, src] = getControl(control, W)
W_out = W;
[bc, src] = deal([]);

% also allow passing control as mat-file
if ischar(control)
    assert(isfile(control));
    tmp = load(control);
    control = tmp.control;
end

if ~isempty(control)
    if isfield(control, 'bc')
        bc = control.bc;
    end
    if isfield(control, 'src')
        src = control.src;
    end
    if isfield(control, 'W')
        % check compatibility
        %assert(numel(W)==numel(control.W) && ...
        %       all(arrayfun(@(x,y)strncmp(x.name, y.name), W, control.W)), ...
        %       'Input well not compatible with realization well');
           
        flds = {'type', 'val', 'sign', 'status', 'cstatus', 'lims'};
        flds = flds(isfield(control.W, flds));
        for k = 1:numel(W)
            for j = 1:numel(flds)
                W(k).(flds{j}) = control.W(k).(flds{j});
            end
        end
    end
end
end
 
%--------------------------------------------------------------------------


function h = getHandlers(outDir, fields)
[dd, fldr] = fileparts(outDir);
for k = 1:numel(fields)
    h.(fields{k}) = ResultHandler('dataDirectory', dd, ...
        'dataFolder', fldr, 'writeToDisk', true, ...
        'dataPrefix', fields{k} );
end
end

%--------------------------------------------------------------------------

function W = adjustControlsForNonUnitBFactors(W, model, p_ref)
% we assume all injectors are injecting water for now, check that this is
% the case
iix      = vertcat(W.sign) > 0;
compi    = vertcat(W.compi);
compiInj = all(compi(iix, :));

if ~(model.water && (size(compiInj, 1)==1 || all(compiInj(2:end)==0)))
    warning('Injection of other fluids than water is not currently handled correctly here, update ...')
end

bW = model.fluid.bW(p_ref); % pressure should be immaterial here
if abs(bW-1) > 1e-5
    for k = 1:numel(W)
        if strcmp(W(k).type, 'rate')
            W(k).val = W(k).val/bW;
        end
    end
end
end

%--------------------------------------------------------------------------

% function state = adjustSolutionForNonUnitBFactors(state, model, p_ref)
% bW = model.fluid.bW(p_ref); % pressure should be immaterial here
% if abs(bW-1) > 1e-5
%     for k = 1:numel(state.wellSol)
%         state.wellSol.flux = state.wellSol.flux*bW;
%     end
% end
% end
