function partials = perturbWell(model, state, W, getResiduals, posControl, varargin)
% Undocumented utility function used for adjoint-based trajectory optimization 

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

opt = struct('approxType',   'twoSided');
opt = merge_options(opt, varargin{:});

% partials: struct containing
%   dFdu      : cell-array containing partial derivative approximations such
%               that dFdu{eqno}(:,k) is partial derivative of equation eqno wrt to parameter i
%   pointIx   : parameter to point index (parameter i -> point pntIx(i))
%   direction : unit perturbation directions (parameter i -> direction dir(i,:))
%
% try calling
r = getResiduals(state, W);
nc = numel(state.pressure);
[pntIx, v] = deal(posControl.parameters.pointIx, posControl.parameters.perturbation);
nu = numel(pntIx); % number of parameters

sysNms = fieldnames(r); % system nammes (pressure, backward, forward)
for k = 1:numel(sysNms)
    if ~isempty(r.(sysNms{k}))
        neq = numel(r.(sysNms{k}).eqs);
        partials.(sysNms{k}) = repmat({cell(1, nu)}, [1, neq]);
    else
        partials.(sysNms{k}) = [];
    end
end

% get index of well
wno  = strcmp(posControl.w.name, {W.name});
assert(nnz(wno)==1, '')


[w, ws] = deal(W(wno), state.wellSol(wno));
pertFac = .5;
if strcmp(opt.approxType, 'oneSided')
    pertFac = 1;
    r_m    = r;
    wcIx_m = sparse(W(wno).cells, 1, true, nc, 1);
end
        
points0 = posControl.controlPoints;
for uno = 1:nu
    pert = v(uno,:);
    %traj_p = posControl.controlPoints;
    %traj_p(pntIx(uno), :) = traj_p(pntIx(uno), :) + pertFac*pert;
    % reset
    posControl.controlPoints = points0;
    posControl.controlPoints(pntIx(uno),:) = points0(pntIx(uno), :) + pertFac*pert;
    %posControl.controlPoints
    %[w_p, ws_p] = posControl.perturbFun(w, ws, traj_p);
    [w_p, ws_p] = updateWellTrajectory(model, w, ws, posControl.getTrajectory());
    if ~w_p.status
        warning('Well-placement outside grid, this is currently not handled');
        wcIx_p = sparse(nc, 1);
    else
        wcIx_p = sparse(w_p.cells, 1, true, nc, 1);
    end
    [W(wno), state.wellSol(wno)] = deal(w_p, ws_p);
    r_p = getResiduals(state, W); % residuals of perturbed eqs pos
    %r_p = rr.(nm).eqs;
    if strcmp(opt.approxType, 'twoSided')
        %traj_m = posControl.controlPoints;
        %traj_m(pntIx(uno), :) = points0(pntIx(uno), :) - pertFac*pert;
        posControl.controlPoints = points0;
        posControl.controlPoints(pntIx(uno),:) = points0(pntIx(uno), :) - pertFac*pert;
        %posControl.controlPoints
        [w_m, ws_m] = updateWellTrajectory(model, w, ws, posControl.getTrajectory());
        if ~w_m.status
            wcIx_m = sparse(nc, 1);
        else
            %[w_m, ws_m] = posControl.perturbFun(w, ws, traj_m);
            wcIx_m = sparse(w_m.cells, 1, true, nc, 1);
        end
        [W(wno), state.wellSol(wno)] = deal(w_m, ws_m);
        r_m = getResiduals(state, W); % residuals of perturbed eqs neg
        %r_m = rr.(nm).eqs;
    end
    % loop through systems
    for sno = 1:numel(sysNms)
        sysNm = sysNms{sno};
        if ~isempty(r.(sysNm))
            % loop through system equations
            [eqs, eqTps] = deal(r.(sysNm).eqs, r.(sysNm).eqTypes);
            for eqno = 1:numel(eqs)
                if strcmp(eqTps{eqno}, 'cell')
                    % perturbation should only hav effect on well-cells  
                    cix = find(wcIx_p | wcIx_m);
                    tmp = (r_p.(sysNm).eqs{eqno}(cix)-r_m.(sysNm).eqs{eqno}(cix))/norm(pert);
                    partials.(sysNm){eqno}{uno} = sparse(cix, 1, tmp, nc ,1);
                    %partials.(sysNm){eqno}{uno} = (r_p.(sysNm).eqs{eqno}-r_m.(sysNm).eqs{eqno})/norm(pert);
                else
                    % perturbed might be outside grid -> well becomes inactive
                    r_p_tmp = checkStatus(r_p.(sysNm).eqs{eqno}, W, wno);
                    r_m_tmp = checkStatus(r_m.(sysNm).eqs{eqno}, W, wno);
                    partials.(sysNm){eqno}{uno} = (r_p_tmp-r_m_tmp)/norm(pert);
                end
            end
        end
    end
end
% concatenate over parameters
for sno = 1:numel(sysNms)
    if ~isempty(partials.(sysNms{sno}))
        partials.(sysNms{sno}) = cellfun(@(x)horzcat(x{:}), partials.(sysNms{sno}), 'UniformOutput', false);
    end
end
% reset
posControl.controlPoints = points0;
end

function r = checkStatus(r, W, wno)
if numel(r) ~= numel(W)
    assert(numel(r)+1 == numel(W), 'Check what''s going on here!');
    tmp = zeros(numel(W), 1);
    tmp(~wno(:)) = r;
    r = tmp;
end
end
    
