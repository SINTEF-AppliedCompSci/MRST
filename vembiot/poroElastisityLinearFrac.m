function states = poroElastisityLinearFrac(state0, G, problem, schedule, varargin)%bc_f, bc_s)
%function states = poroElastisityLinear(G, disc, rock, fluid, schedule, varargin)%bc_f, bc_s)

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
%}
opt = struct('do_plot', true, ...
    'linsolve', @mldivide, ...
    'disc', 'tpfa','tens_strengt',1e8,'do_patch',true,'solve_meth','implicit',...
    'min_E',1e6,...
    'max_perm',1*darcy,...
    'alpha',0,...
    'nobiot',false,...
    'max_fcells',1e8);
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
well_cells=problem.rock.perm(:,1)>1e-2*darcy;%(G.cells.num,1);
C(well_cells,:)=C(well_cells,:)*1e-6;
[uu_t, extra] = VEM_linElast(G, C, problem.el_bc, problem.load, 'dual_field', false);
%[S,op]=VEM_assemble(G,C,'add_operators',true);
D=extra.D;
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
switch opt.disc
    case 'tpfa'
        T = computeTrans(G, problem.rock);
        state_f = lincompTPFA(dt, state_f, G, T, pv, problem.fluid,problem.rock,'MatrixOutput', true,...
            'wells', problem.W,'src', problem.src , 'bc', problem.bc_f);
    case 'mpfa'
        T = computeMultiPointTrans(G, problem.rock);
        state_f = lincompMPFA(dt, state_f, G, T, pv, problem.fluid, 'MatrixOutput', true,...
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
frac_cells=false(size(well_cells));%(G.cells.num,1);
fbc = addFluidContribMechVEM(G, problem.bc_f, problem.rock, isdirdofs);
for i = 1:numel(schedule.step.val)
    dt = schedule.step.val(i);
    not_converged=true;
    org_frac_cells=frac_cells;
    count_it=1;max_it=10;
    t = t+dt;
    %rhs_s do not change
    x_prev = x;
    p_prev=p;
    while not_converged && (count_it < max_it)     
     
        mat = sparse(size(Af, 1)-G.cells.num, size(div, 2));
        fac = sparse(1:G.cells.num, 1:G.cells.num, problem.rock.alpha, G.cells.num, G.cells.num);
        fac2 =fac;
        matdp=[-gradP*fac,mat'];
        matpd=[fac*div;mat];
        switch opt.solve_meth
            case 'seq_pressure'
                %fac = sparse(1:G.cells.num, 1:G.cells.num, problem.rock.alpha, G.cells.num, G.cells.num);
                rhsf(1:G.cells.num) = orhsf(1:G.cells.num)*dt+ct(1:G.cells.num, 1:G.cells.num)*p_prev;
                rhsf(G.cells.num+1:end) = orhsf(G.cells.num+1:end);
                xp_g=opt.linsolve(ct+dt*Af,rhsf);
                xd_prev=x(1:numel(rhs_s));
                xd=opt.linsolve(As,rhs_s-fbc-matdp*xp_g);
                %xp=opt.linsolve(ct+dt*Af,rhsf-matpd*(xd-xd_prev)); % do
                %work
                %not 
                xp=xp_g;
                x=[xd;xp]; 
            case 'seq_stress'
                %fac = sparse(1:G.cells.num, 1:G.cells.num, problem.rock.alpha, G.cells.num, G.cells.num);
                rhsf(1:G.cells.num) = orhsf(1:G.cells.num)*dt+ct(1:G.cells.num, 1:G.cells.num)*p_prev;
                rhsf(G.cells.num+1:end) = orhsf(G.cells.num+1:end);
                %xp=opt.linsolve(ct+dt*Af,rhsf);
                xp_prev=x(numel(rhs_s)+1:end);
                xd_prev=x(1:numel(rhs_s));
                xd=opt.linsolve(As,rhs_s-fbc-matdp*xp_prev);
                xp=opt.linsolve(ct+dt*Af,rhsf-matpd*(xd-xd_prev));
                x=[xd;xp];
                % do not work ??
            case 'imp_frac'
                smod=0;
                rhsf(1:G.cells.num) = orhsf(1:G.cells.num)*dt+ct(1:G.cells.num, 1:G.cells.num)*p_prev+smod*fac*div*x_prev(ind_s);% can be moved out
                rhsf(G.cells.num+1:end) = orhsf(G.cells.num+1:end);
                rhs = [rhs_s-fbc;
                    rhsf];
                SS = [As, matdp;...
                    smod*matpd, ct+dt*Af];
                x = opt.linsolve(SS, rhs);
                %{
                ind_frac=false(size(x));
                ind_frac(ind_s)=true;
                ind_frac(ind_f(frac_cells))=true;
                %}                
                %error()
            case 'implicit'
                %fac = sparse(1:G.cells.num, 1:G.cells.num, problem.rock.alpha, G.cells.num, G.cells.num);
                if(opt.nobiot)
                  fac2 = sparse(1:G.cells.num, 1:G.cells.num, problem.rock.alpha.*frac_cells, G.cells.num, G.cells.num);
                  matpd=[fac2*div;mat];
                end
                    
                
                rhsf(1:G.cells.num) = orhsf(1:G.cells.num)*dt+ct(1:G.cells.num, 1:G.cells.num)*p_prev+fac2*div*x_prev(ind_s);
                rhsf(G.cells.num+1:end) = orhsf(G.cells.num+1:end);
                rhs = [rhs_s-fbc;
                    rhsf];
                %{
                SS = [As, [(-gradP)*fac, mat'];...
                    [fac*div; mat], ct+dt*Af];
                    %}
                SS = [As, matdp;...
                    matpd, ct+dt*Af];    
                %a = tic;
                x = opt.linsolve(SS, rhs);
            otherwise
                error()
                
        end
        %a = toc(a);
        
        p = x(ind_f);
        u(isdirdofs) = Vdir(isdirdofs);
        u(~isdirdofs) = x(ind_s);
        uu = reshape(u, G.griddim, [])';
        %%{
        % fracture calc
        if(G.griddim==2)
            nlin=3;
        else
            nlin=6;
        end
        % from assambly
        %[S,op]=VEM_assemble(G,C,'add_operators',true);
       
        stress=reshape(D*extra.WC'*extra.assemb'*reshape(uu',[],1),3,[])';
        %strain=reshape(op.WC'*op.assemb'*reshape(uu',[],1),3,[])';
        stress=calculateStressVEM(G,uu, extra);
        if(opt.do_patch)
            stress=patchRecovery(G,stress,'filter',~frac_cells & ~well_cells);  
        end
        [sigm,evec]=calStressEigsVEM(G,stress);
        %
        %cells=find(sigm(:,2)>opt.tens_strengt & ~frac_cells);
        cells=find(sigm(:,2)>opt.tens_strengt & ~frac_cells & ~well_cells);
        if(numel(cells)>0)
            %[mm,i]=min(sigm(cells,1));cells=cells(i);
            [d,i]=sort(sigm(cells,2),'descend');
            if(opt.max_fcells<numel(cells))
                cells=cells(i(1:opt.max_fcells));
            else
                cells=cells(i);
            end
            %
            nc=numel(cells);
            C_tmp=Enu2C(opt.min_E*ones(nc,1),0.3*ones(nc,1),struct('griddim',G.griddim));
            C(cells,:)=C_tmp;
            
            D_all = C2D(C, G);
            mcells=1:G.cells.num;
            D = reshape(D_all(mcells, :)', nlin, [])';
            [i, j] = blockDiagIndex(nlin*ones(numel(mcells), 1), nlin*ones(numel(mcells), 1));
            D = sparse(i, j, reshape(D', [], 1), nlin*numel(mcells), nlin*numel(mcells));
            
            
            frac_cells(cells)=true;
            %[~, extra] = VEM_linElast(G, C, problem.el_bc, problem.load, 'dual_field', false);
            %As = extra.disc.A;
            %S=extra.assemb*op.WC*D*op.WC'*extra.assemb';
            KH = extra.volmap*extra.WC*D*extra.WC' + (extra.I - extra.PP)'*extra.SE*(extra.I - extra.PP);
            S=extra.assemb*KH*extra.assemb';
            % correct if direclet condtions is zero
            As   = S(~extra.disc.isdirdofs, ~extra.disc.isdirdofs);
            %
            problem.rock.perm(cells,:)=opt.max_perm*ones(nc,1);
            T = computeTrans(G, problem.rock);
            state_f = lincompTPFA(dt, state_f, G, T, pv, problem.fluid,problem.rock, 'MatrixOutput', true,...
                'wells', problem.W,'src', problem.src , 'bc', problem.bc_f);
            Af=state_f.A;
        end
        figure(99)
        clf,%plotGrid(G),plotGrid(G,find(frac_cells),'FaceColor','r')
        val=sigm(:,2);
        val(frac_cells)=min(val(frac_cells));
        plotCellData(G,val)
        plotGrid(G,frac_cells,'FaceColor','k');colorbar
        pause(0.01)
        %disp(min(sigm(:,1)));
        %{
        figure(99)
        clf,%plotGrid(G),plotGrid(G,find(frac_cells),'FaceColor','r')
        val=sigm(:,1);
        val(frac_cells)=min(val(frac_cells));
        plotCellData(G,val)
        plotGrid(G,frac_cells,'FaceColor','k');
        figure(97)
        clf,plotCellData(G,p)
        figure(98)
        clf,plotCellDataDeformed(G,sigm(:,1),uu),colorbar
        %pause()
        %}
        not_converged=not(all(org_frac_cells == frac_cells));
        if(not_converged)
            %x=x_prev;
            count_it=count_it+1;
        end
    end
    states{count} = struct('t', t, 'pressure', p, 'uu', uu, 'bhp', x(ind_f(end)+1:end),'frac_cells',frac_cells,'count_it',count_it);
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




