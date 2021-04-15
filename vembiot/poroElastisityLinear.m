function states = poroElastisityLinear(state0, G, problem, schedule, varargin)%bc_f, bc_s)
%function states = poroElastisityLinear(G, disc, rock, fluid, schedule, varargin)%bc_f, bc_s)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
%}
opt = struct('do_plot', true, ...
    'linsolve', @mldivide, ...
    'disc', 'tpfa','eta',0);
%'type', 'mrst', 'linsolve', @mldivide, ...
%    'local_assembly', false, 'do_local', falsep
%    'dual_field', true);
opt = merge_options(opt, varargin{:});
%assert(G.griddim == 2)
%el_bc.displacement = testcase.disp_sol(G.nodes.coords);
%[uu_t, p_t, S_t, A, extra] = VEM2D_linElast2f(G, problem.Ev, problem.nuv,
%problem.el_bc, problem.load, 'dual_field', false);
%C = Enu2C(problem.Ev, problem.nuv,G);
C = problem.C;
[uu_t, extra] = VEM_linElast(G, C, problem.el_bc, problem.load, 'dual_field', false);
As = extra.disc.A;
gradP = extra.disc.gradP;
div = extra.disc.div;
isdirdofs = extra.disc.isdirdofs;
rhs_s = extra.disc.rhs;
Vdir = extra.disc.V_dir;
ind_s = [1:size(As, 1)]';%#ok

% calculate discretization of flow

pv = poreVolume(G, problem.rock);
pressure = 100*barsa*ones(G.cells.num, 1);
state_f = struct('pressure', pressure, 's', ones(G.cells.num, 1), 'flux', zeros(G.faces.num, 1));
dt = schedule.step.val(1);%fake
fluid_tmp=problem.fluid;
switch opt.disc
    case 'tpfa'
        T = computeTrans(G, problem.rock);
        state_f = lincompTPFA(dt, state_f, G, T, pv, fluid_tmp, 'MatrixOutput', true,...
            'wells', problem.W,'src', problem.src , 'bc', problem.bc_f);
    case 'mpfa'
        T = computeMultiPointTrans(G, problem.rock,'eta',opt.eta);
        state_f = lincompMPFA(dt, state_f, G, T, pv, fluid_tmp, 'MatrixOutput', true,...
            'wells', problem.W,'src', problem.src , 'bc', problem.bc_f);
        warning('some boundaries may be wrong')
    otherwise 
       error('No such discretization implemented')
end

        % get discretization
% definitions with out dt
Af = state_f.A;
orhsf = state_f.rhs;
ct = state_f.ct;
ind_f = [ind_s(end)+1:ind_s(end)+G.cells.num]';%#ok

% initial variables
x = zeros(ind_f(end), 1);% set all to zero first
x(ind_f) = state0.pressure;% set cell pressure part
p = state0.pressure; % variable for pressure
u_tmp = reshape(state0.uu', [], 1);% initial displacement
x(1:ind_s(end)) = u_tmp(~isdirdofs);
u = u_tmp;


rhsf = zeros(size(orhsf));
plotops = {'EdgeColor', 'none'};
count = 1;
t = 0;
states = cell(numel(schedule.step.val), 1);
for i = 1:numel(schedule.step.val)
    dt = schedule.step.val(i);
    t = t+dt;
    %rhs_s do not change
    xo = x;
    fac = sparse(1:G.cells.num, 1:G.cells.num, problem.rock.alpha, G.cells.num, G.cells.num); 
    rhsf(1:G.cells.num) = orhsf(1:G.cells.num)*dt+ct(1:G.cells.num, 1:G.cells.num)*p+fac*div*x(ind_s);
    rhsf(G.cells.num+1:end) = orhsf(G.cells.num+1:end);
    % fake boundary 
    
    
    % set boundary pressure to 100 barsa need to be taken more care full
    % need to consider best way of handaling boundary condtion to get
    % symetric system ???, pressure calculations on boundary seems to be
    % nessesary to calculate correct gradient for pressure contribution to
    % force for the elasiticyt part
    %NB is use of fac correct ???
    fbc = addFluidContribMechVEM(G, problem.bc_f, problem.rock, isdirdofs);
    
    
    
    %assert(max(max(abs(sum(gradP*fac, 2))))<eps*1e3);
    % can use the fix from the ad solvers addBcFromFluid(model, bc)

    %rhs = [rhs_s-sum(gradP*fac, 2)*p_bc;
    %    rhsf];
    rhs = [rhs_s-fbc;
        rhsf];
    % to add bhp wells as separate variables
    %fac = 0;
    mat = sparse(size(Af, 1)-G.cells.num, size(div, 2));
    SS = [As, [(-gradP)*fac, mat'];...
        [fac*div; mat], ct+dt*Af];
    %a = tic;
    x = opt.linsolve(SS, rhs);
    %a = toc(a);
    
    p = x(ind_f);
    u(isdirdofs) = Vdir(isdirdofs);
    u(~isdirdofs) = x(ind_s);
    uu = reshape(u, G.griddim, [])';   
    states{count} = struct('t', t, 'pressure', p, 'uu', uu, 'bhp', x(ind_f(end)+1:end));
    count = count+1;
    if(opt.do_plot)
        figure(1), clf
        subplot(2, 2, 1), cla
        plotNodeData(G, uu(:, 1), plotops{:}), colorbar;
        subplot(2, 2, 2), cla
        plotNodeData(G, uu(:, 2), plotops{:});colorbar
        subplot(2, 2, 3), cla
        plotCellDataDeformed(G, p/barsa, uu);colorbar
        subplot(2, 2, 4), cla
        ovdiv = extra.disc.ovol_div;
        mdiv = ovdiv*reshape(uu', [], 1)./G.cells.volumes;
        plotCellDataDeformed(G, mdiv, uu);colorbar(),
        %
        figure(2), clf
        uur = uu-state0.uu;
        subplot(2, 2, 1), cla
        plotNodeData(G, uur(:, 1), plotops{:}), colorbar;
        subplot(2, 2, 2), cla
        plotNodeData(G, uur(:, 2), plotops{:});colorbar
        subplot(2, 2, 3), cla
        plotCellDataDeformed(G, p/barsa, uur);colorbar
        subplot(2, 2, 4), cla
        ovdiv = extra.disc.ovol_div;
        mdiv = ovdiv*reshape(uur', [], 1)./G.cells.volumes;
        plotCellDataDeformed(G, mdiv, uu);colorbar(),
        pause(0.01);
    end
end
end




