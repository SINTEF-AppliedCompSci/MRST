classdef DiagnosticsObjective < handle
    % Base class for diagnostics objectives, does not contain a particular
    % objetive
    properties
        model                          % diagnostics models
        wellControl       = true;      % compute derivatives wrt well controls
    end
    
    methods
        function obj = DiagnosticsObjective(model, varargin)
            obj = merge_options(obj, varargin{:});
            obj.model = model; 
            obj.model.parentModel.parentModel.G = addBoundingBoxFields(obj.model.parentModel.parentModel.G);
        end
        % -----------------------------------------------------------------
        
        function result = compute(obj, W, varargin)
            opt = struct('state0', obj.model.state0, ...
                         'state',              [], ...  % if already computed
                          'D',                 [], ...  %
                          'computeGradient',   true, ...
                          'outputState',       true, ...
                          'outputDiagnostics', true);
            opt = merge_options(opt, varargin{:});
            result = struct('value', [], 'gradient', [], 'state', [], 'D', []);
            [state0, state, D] = deal(opt.state0, opt.state, opt.D);
            
            forces = obj.model.getValidDrivingForces();
            forces.W = W;
            obj.model  = obj.model.validateModel(forces);
            
            % solve pressure equation
            evalDiagn = false;
            if isempty(state)
                 state = obj.model.solvePressure(W, 'state0', state0);
                 evalDiagn = true;
            end
            % solve diagnsotics
            if isempty(D) || evalDiagn
                 [state, D] = obj.model.solveDiagnostics(state);
            end
            if opt.outputState
                result.state = state;
            end
            if opt.outputDiagnostics
                result.D = D;
            end
            % Objective J
            %v = obj.evaluate(state0, state, D, computeGradient);
                
            % Adjoint:  dFdx'*L = -dJdx
            % Gradient: dFdu'*L + dJdu
            if ~opt.computeGradient
                result.value = obj.evaluate(state, D, W);
            else
                [result.value, dJdx, dJdu] = obj.evaluate(state, D, W);
                % solve adjoint system
                system = obj.model.getCoupledSystem(state, forces);
                L = obj.model.solveAdjoint(system, dJdx);
                
                %%%
                L.backward = {L.backward.tof};
                %%%
                % get gradient
                dFdu = obj.getEquationPartials(system, state, forces);
                % sum up contributions
                g = struct('well', [], 'position', {cell(1, numel(W))});
                if obj.wellControl
                    g.well = assembleGradient(L, dFdu.well, dJdu.well, W, true);
                end
                if isfield(W(1), 'posControl')
                    pcontrols = {W.posControl};
                    posWells = find(~cellfun(@isempty, pcontrols));
                    for k = 1:numel(posWells)
                        ci = posWells(k);
                        gk = assembleGradient(L, dFdu.position{ci}, dJdu.position{ci}, W, false);
                        %g.position{ci} = assembleGradient(L, dFdu.position, dJdu.position);
                        % map gradient to 3D
                        
                        pc = pcontrols{ci};
                        v  = pc.parameters.perturbation;
                        % normalize
                        v  = bsxfun(@rdivide, v, sqrt(sum(v.^2,2)));
                        % expand
                        gk = bsxfun(@times, gk, v);
                        % sum over pointIx
                        np = max(pc.parameters.pointIx);
                        tmp = cell(np, 1);
                        for pk = 1:np
                            tmp{pk} = sum( gk(pc.parameters.pointIx == pk,:) ).';
                        end
                        g.position{ci} = vertcat(tmp{:}); 
                    end
                end
                result.gradient = g;
            end
        end
        % -----------------------------------------------------------------
        
        function varargout = evaluate(obj, state, D, W)  %#ok
            error('''DiagnosticsObjective'' does not contain an objective that can be evaluated')
        end
        % -----------------------------------------------------------------
        
        function dFdu = getEquationPartials(obj, system, state, forces)
            dFdu = struct('well', [], 'position', []);
            obj.model = obj.model.validateModel(forces);
            if obj.wellControl
                % only non-zero for closureEq
                neq = numel(system.pressure.eqs);
                dFdu.well.pressure = repmat({[]}, [1, neq]);
                ix = strcmp(system.pressure.eqNames, 'closureWells');
                assert(nnz(ix)==1);
                sc = getControlEquationScalings(system.pressure, forces.W);
                dFdu.well.pressure{ix} = -spdiags(sc, 0, numel(sc), numel(sc));
                
                %dFdu.well.pressure{ix} = -speye(nnz([forces.W.status]));
                [dFdu.well.forward, dFdu.well.backward] = deal([]);
            end
            W = forces.W;
            if isfield(W(1), 'posControl')
                pcontrol = {W.posControl};
                wellNo = find(~cellfun(@isempty, pcontrol));
                if any(wellNo)
                    % get reference psystem (also available in system ...)
                    psys = obj.model.getPressureSystem(state, forces, true);
                    [neq, nc] = deal(numel(psys.eqs), obj.model.G.cells.num);
                    getResiduals = @(state, W)getSystemResiduals(obj.model, state, W);
                    %eqTypes = system.pressure.eqTypes;
                    for k = 1:numel(wellNo)
                        cno = wellNo(k);
                        %for cno = 1:numel(obj.positionControl)
                        gboModel = obj.model.parentModel.parentModel;
                        %if forces.W(cno).status
                        dFdu.position{cno} = perturbWell(gboModel, state, forces.W, getResiduals, pcontrol{cno});                        
                    end
                end
            end
        end
    end
end
% -----------------------------------------------------------------
% -----------------------------------------------------------------

function r = getSystemResiduals(model, state, W)
% this can be worked on to get more efficient, e.g., modify model directly
forces = model.getValidDrivingForces();
forces.W = W;
r = model.getCoupledSystem(state, forces, 'resOnly', true, 'validate', false);
end
% -----------------------------------------------------------------

function g = assembleGradient(L, dFdu, dJdu, W, isWellControl)
% check that L and dFdu contain same non-empty flds
flds = {'pressure', 'forward', 'backward'};
g = 0;
for k = 1:numel(flds)
    [Li, Pi] = deal(L.(flds{k}), dFdu.(flds{k}));
    if ~isempty(Pi)
        for eqno = 1:numel(Pi)
            if ~isempty(Pi{eqno})
%                 if strcmp(flds{k}, 'pressure') && ~strcmp(system.(flds{k}).eqTypes{eqno}, 'cell')
%                     tmp = double(stat);
%                     tmp(stat) = Li{eqno};
%                     g = g + Pi{eqno}'*tmp;
%                 else
                    g = g + Pi{eqno}'*Li{eqno};
                %end
            end
        end
    end
end
if ~isempty(dJdu)
    g = g + vertcat(dJdu);
end
stat = vertcat(W.status);
if ~all(stat) && isWellControl
    tmp = double(stat);
    tmp(stat) = g;
    g = tmp;
end
end
% -----------------------------------------------------------------

function s = getControlEquationScalings(sys, W)
s = ones(nnz([W.status]),1);
% account for scaled bhp-control-eqs
bhpcnt = strcmp('bhp', {W.type});
if any(bhpcnt)
    ecno = strcmp('closureWells', sys.eqNames);
    ec   = sys.eqs{ecno};
    assert(isa(ec, 'GenericAD'));
    vno  = find(strcmp('bhp', sys.varNames));
    jno  = find(ec.offsets <= vno, 1, 'last');
    J    = ec.jac{jno};
    if isa(J, 'DiagonalJacobian')
        tmp = J.diagonal(:, vno-ec.offsets(jno)+1);
    else
        tmp = spdiags(J);
    end
end
s(bhpcnt) = tmp(bhpcnt);
end


function G = getGridWithBoundingBoxes(obj)
G = obj.model.parentModel.parentModel.G;
if ~isfield(G.faces, 'bbox')    
    G = computeBoundingBoxes(G);
end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
