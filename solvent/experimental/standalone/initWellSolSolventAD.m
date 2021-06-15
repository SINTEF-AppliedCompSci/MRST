function wellSol = initWellSolAD(W, model, state0, wellSolInit)
%Set up well solution struct for a automatic differentiation model
%
% SYNOPSIS:
%   wellSol = initWellSolAD(W, model, state0);
%   wellSol = initWellSolAD(W, model, state0, ws);
%
% DESCRIPTION:
%   Create or extract the wellSol, and ensure that it contains the correct
%   fields for advanced solvers with well limits and variable perforation
%   counts. This function will first look for a explicitly passed wellSol
%   to modify, then it will consider any wellSol residing in state0. If
%   neither is found, it will attempt to construct one based on W.
%
% REQUIRED PARAMETERS:
%   W          - Control well for which we are going to create a well
%                solution structure.
%
%   model      - Subclass of ReservoirModel. Used to determine how many
%                and which phases are present.
%
%   state0     - State, possibly with a wellSol given already (see
%                initResSol/initState).
%
%   wellSolInit - Initial wellSol.
%
%
% RETURNS:
%   wellSol     - Well solution struct with additional fields ready for
%                 simulation.
%

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

wellSolGiven = (nargin == 4);

if wellSolGiven
    wellSol = wellSolInit;
elseif isfield(state0, 'wellSol')
    wellSol = state0.wellSol;
else
    wellSol = defaultWellSol(state0, W, model);
end
wellSol = assignFromSchedule(W, model, wellSol);
end

function ws = defaultWellSol(state, W, model)
nw = numel(W);
[~, actPh] = model.getActivePhases();

ws = repmat(struct(...
    'name',   [],...
    'status', [],...
    'type',   [],...
    'val',    [],...
    'sign',   [],...
    'bhp',    [],...
    'qWs',    [],...
    'qOs',    [],...
    'qGs',    [],...
    'mixs',   [],...
    'cstatus',[],...
    'cdp',    [],...
    'cqs',    []), [1, nw]);

 
% Additional model dependent fields
if isprop(model, 'solvent') && model.solvent % solvent model
	 [ws(:).qSs] = deal(0);
end

if isprop(model, 'polymer') && model.polymer % polymer model
	 [ws(:).qWPoly] = deal(0);
end

if isprop(model, 'surfactant') && model.surfactant % surfactant model
	 [ws(:).qWSft] = deal(0);
end

if isprop(model, 'compositionalFluid') % Compositional model
     ncomp = model.compositionalFluid.getNumberOfComponents();
	 [ws(:).components] = deal(zeros(1, ncomp));
end

% just initialize fields that are not assigned in assignFromSchedule
for k = 1:nw
    nConn = numel(W(k).cells);
    nPh   = numel(actPh);
    ws(k).name = W(k).name;
    % To avoid switching off wells, we need to start with a bhp that makes
    % a producer produce and an injector inject. Hence, we intitialize the
    % bhp such that the top connection pressure is 5bar above/below the
    % corresponding well-cell pressure. If W(k).dZ is ~= 0, however, we
    % don't know wht a decent pressure is ...
    % The increment should depend on the problem and the 5bar could be a
    % pit-fall... (also used in initializeBHP in updateConnDP)
    ws(k).bhp = state.pressure(W(k).cells(1)) + 5*W(k).sign*barsa;

    irate = eps;
    if model.water
        ws(k).qWs  = W(k).sign*irate;
    end
    if model.oil
        ws(k).qOs  = W(k).sign*irate;
    end
    if model.gas
        ws(k).qGs  = W(k).sign*irate;
    end
    if isprop(model, 'solvent') && model.solvent
       ws(k).qSs = W(k).sign*irate;
    end
    if isprop(model, 'polymer') && model.polymer
       ws(k).qWPoly = W(k).poly*ws(k).qWs;
    end
    if isprop(model, 'surfactant') && model.surfactant
       ws(k).qWSft = W(k).surfact*ws(k).qWs;
    end
    
    ws(k).mixs = W(k).compi;
    ws(k).qs   = W(k).sign*ones(1, nPh)*irate;
    ws(k).cdp  = zeros(nConn,1);
    ws(k).cqs  = zeros(nConn,nPh);
end
end

function ws = assignFromSchedule(W, model, ws)
% set fields that should be updated if control has changed
for k = 1:numel(W)
    ws(k).status  = W(k).status;
    ws(k).type    = W(k).type;
    ws(k).val     = W(k).val;
    ws(k).sign    = W(k).sign;
    ws(k).cstatus = W(k).cstatus;

    tp = W(k).type;
    v  = W(k).val;

    switch tp
        case 'bhp'
            ws(k).bhp = v;
        case 'rate'
            if model.water
                ws(k).qWs = v*W(k).compi(1);
            end
            if model.oil
                ix = 1 + model.water;
                ws(k).qOs = v*W(k).compi(ix);
            end
            if model.gas
                ix = 1 + model.water + model.oil;
                ws(k).qGs = v*W(k).compi(ix);
            end
            if isprop(model, 'solvent') && model.solvent
                ix = 1 + model.water + model.oil + model.gas;
                ws(k).qSs = v*W(k).compi(ix);
            end
            if isprop(model, 'polymer') && model.polymer
                ws(k).qWPoly = W(k).poly*ws(k).qWs;
            end
            if isprop(model, 'surfactant') && model.surfactant
               ws(k).qWSft = W(k).surfact*ws(k).qWs;
            end
        case 'orat'
            ws(k).qOs = v;
        case 'wrat'
            ws(k).qWs = v;
        case 'grat'
            ws(k).qGs = v;
        case 'srat'
            ws(k).qSs = v;
    end % No good guess for qOs, etc...
end
end

