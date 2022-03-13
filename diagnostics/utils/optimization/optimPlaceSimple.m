function u = optimPlaceSimple(p, problem, varargin)
% simple gradient-based search for well-trajectory optimization. See 
% e.g., optimizeDiagnosticsPositionEgg.m for usage 

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

opt = struct('stepInit',         1, ...
             'maxSteps',        20, ...
             'maxRelative',       .05, ...
             'maxLineSearchIts', 5, ...
             'plotTrajectories', true, ...
             'plotAx',             []);
% if isempty(p.controlList)
%     W  = p.initialGuess;         
%     u0 = p.convertInputToControl(W);
% else
%     u0 = p.controlList{end};
% end

u0 =p.scaleVariables(p.fetchVariablesFun(problem));

[v, g] = p.getScaledObjective(u0);
u = u0;
cix = (1:numel(u))';
%cix = (1:p.n.well)';
%pix = p.n.well+ (1:p.n.pos)';
fprintf('Initial objective value: %e\n', v)
% get actual u used
%u = p.convertInputToControl(p.initialGuess);
step = opt.stepInit;

pc = {p.W.posControl};
pc = horzcat(pc{:});
if opt.plotTrajectories 
    if isempty(opt.plotAx)
        figure, hold on
        ax = gca;
        plotGrid(pc(1).G, 'FaceColor', 'none', 'EdgeColor', [.6 .6 .6], 'EdgeAlpha', .4);
        view(3), drawnow;
    else
        ax = opt.plotAx;
    end
end
%optIx = 1;     
%pc = pc(1)
for k =1:opt.maxSteps
    fprintf('Outer iteration: %d\n', k);
    ok = false;
    lits = 0;
    if opt.plotTrajectories
        for kk = 1:numel(pc)
            traj = pc(kk).getTrajectory;
            plot3(ax, traj(:,1), traj(:,2), traj(:,3),'-g','LineWidth', 2);
            pnts = pc(kk).controlPoints;
            plot3(ax, pnts(:,1), pnts(:,2), pnts(:,3),'ok','LineWidth', 2);
        end
        drawnow
    end
    
    while ~ok && lits < opt.maxLineSearchIts
        du_cur = step*g;
        % clip according to position constraints
        du_cur = (2^(-lits))*getProjectedUpdate(pc, u,du_cur);
        % clip maxRealitve for controls
        du_c   = sign(du_cur(cix)).*min(opt.maxRelative, abs(du_cur(cix)));
        u_c    = max(0, min(1, u(cix) + du_c(cix)));
        du_c   = u_c - u(cix);
        du_cur(cix) = du_c;

        % compute new objective
        [v_test, g_test] = p.getScaledObjective(u+du_cur);
        fprintf('Computed objective value: %e\n', v_test)
        if opt.plotTrajectories
            for kk = 1:numel(pc)
                traj = pc(kk).getTrajectory;
                plot3(traj(:,1), traj(:,2), traj(:,3),'-r','LineWidth', 2);
            end
            drawnow
        end
        if v_test > v
            ok    = true;
            optIx = numel(p.iterationControls.getValidIds);
        else
           %c = pc.control2param(u);
        end
        lits = lits+1;
    end
    
    if lits == opt.maxLineSearchIts && ~ok
        fprintf('Reached maximal number of line search iterations, exiting.\n')
        break;
    else
        v = v_test;
        g = g_test;
        step = max(.01, norm(du_cur)/norm(g));
        u = u + du_cur;
    end
end
u = p.unscaleVariables(u);
end


function dup = getProjectedUpdate(pc, u, du)
%W = p.convertControlToInput(u);
u_tmp = u+du;
u_tmp = max(0, min(1, u_tmp));
dup = u_tmp - u;
%if p.n.pos > 0
    %pcix = arrayfun(@(w)~isempty(w.posControl), W);
    %pc = {W(pcix).posControl};
    ix = 0;
    for k = 1:numel(pc)
        nParam = pc(k).parameters.nParam;
        ixw = ix + (1:nParam);
        duw = pc(k).getProjectedUpdate(u(ixw), du(ixw), true);
        dup(ixw) = duw;
        ix = ix + nParam;
    end
%end
end