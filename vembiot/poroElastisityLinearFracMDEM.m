function states = poroElastisityLinearFracMDEM(state0, G, problem, schedule, varargin)%bc_f, bc_s)
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
    'max_fcells',1e8,...
    'frac_method','mdem',...
    'e_frac',0.01,...
    'pstep',10);
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
[uu_t, extra] = VEM_linElast(G, C, problem.el_bc, problem.load, 'dual_field', false);%#ok
D=extra.D;
%[S,op]=VEM_assemble(G,C,'add_operators',true);
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
        state_f = lincompTPFA(dt, state_f, G, T, pv, problem.fluid, problem.rock, 'MatrixOutput', true,...
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
%frac_cells=problem.rock.perm(:,1)>1e-2*darcy;%(G.cells.num,1);
e_broke=false(size(G.cells.faces,1),1);% only used when frac_method is mdem+
frac_cells=false(G.cells.num,1);
fbc = addFluidContribMechVEM(G, problem.bc_f, problem.rock, isdirdofs);
for i = 1:numel(schedule.step.val)
    dt = schedule.step.val(i);
    not_converged=true;
    org_frac_cells=frac_cells;
    org_e_broke=e_broke;
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
                xd_prev=x(1:numel(rhs_s));%#ok
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
                error('No such splitting')
                
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
        %{
        mcells=1:G.cells.num;
        D_all = C2D(C, G);
        D = reshape(D_all(mcells, :)', nlin, [])';
        [i, j] = blockDiagIndex(nlin*ones(numel(mcells), 1), nlin*ones(numel(mcells), 1));
        D = sparse(i, j, reshape(D', [], 1), nlin*numel(mcells), nlin*numel(mcells));
        %}
        %stress=reshape(D*extra.WC'*extra.assemb'*reshape(uu',[],1),3,[])';
        %strain=reshape(op.WC'*op.assemb'*reshape(uu',[],1),3,[])';
        %if(opt.do_patch)
        %    stress=patchRecovery(G,stress,'filter',~frac_cells & ~well_cells);
        %end
        %[sigm,evec]=calStressEigs(G,stress);%#ok
        stress=calculateStressVEM(G,uu, extra);
        if(opt.do_patch)
            stress=patchRecovery(G,stress,'filter',~frac_cells & ~well_cells);  
        end
        [sigm,evec]=calStressEigsVEM(G,stress);
        
        %
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
            %C_tmp=Enu2C(opt.min_E*ones(nc,1),0.3*ones(nc,1),struct('griddim',G.griddim));
            %C(cells,:)=C_tmp;
            frac_cells(cells)=true;
            % modify the stiff ness matrix in fracture cells
            [C_tmp,op] = computeDEMC(G,C,'use_dem',true,'cells',find(frac_cells));
            node2dofs =@(nodes) [mcolon(G.griddim*(nodes-1)+1,G.griddim*(nodes-1)+G.griddim)]';
            fndofs=node2dofs(op.nodemap);
            rlen=op.u2rlen*u(fndofs);
            if(G.griddim==2)
                ihf=mcolon(G.cells.facePos(op.cellmap),G.cells.facePos(op.cellmap+1)-1);
                hfn=G.cells.faces(ihf,1);
            else
                error('only valid for 2D sofar')
            end
            %assume simplex
            hfshift=-0.1*(2*(G.faces.neighbors(hfn,1)==rldecode(op.cellmap,(G.griddim+1)*ones(size(op.cellmap))))-1);
            switch opt.frac_method
                case 'mdem'
                    dk_mdem=max(diag(op.k_mdem).*double(rlen<0),opt.min_E);
                    %D_dem=op.Mbl'*diag(dk_mdem)*op.Mbl
                    %DD=reshape(D_dem(op.ind),op.lindim^2,[])';    
                    %C_tmp=C2D(full(DD),G,'inv',true);
                    %C(frac_cells,:)=C_tmp;
                case 'mdem_pluss'
                    e_broke(ihf)=e_broke(ihf) | rlen > opt.e_frac*G.faces.areas(hfn);% break half edge if extension is sufficent.
                    dk_mdem=max(diag(op.k_mdem).*(~(double(rlen>0) & e_broke(ihf))),opt.min_E);
                    
                case 'element'
                    dk_mdem=diag(op.k_mdem)*1e-6;
                otherwise 
                    error('frac_method not implemented')
            end
            D_dem=op.Mbl'*diag(dk_mdem)*op.Mbl;
            DD=reshape(D_dem(op.ind),op.lindim^2,[])';    
            C_tmp=C2D(full(DD),G,'inv',true);
            C(frac_cells,:)=C_tmp;
            %%{
            % recalculated C
            mcells=1:G.cells.num;
            D_all = C2D(C, G);
            D = reshape(D_all(mcells, :)', nlin, [])';
            [i, j] = blockDiagIndex(nlin*ones(numel(mcells), 1), nlin*ones(numel(mcells), 1));
            D = sparse(i, j, reshape(D', [], 1), nlin*numel(mcells), nlin*numel(mcells));
            %}
            %dk_mdem=diag(op.k_mdem)*1e-6;
            %dk_mdem=opt.min_E*ones(size(diag(op.k_mdem)));                
            %[~, extra] = VEM_linElast(G, C, problem.el_bc, problem.load, 'dual_field', false);
            %As = extra.disc.A;
            %S=extra.assemb*op.WC*D*op.WC'*extra.assemb';
            KH = extra.volmap*extra.WC*D*extra.WC' + (extra.I - extra.PP)'*extra.SE*(extra.I - extra.PP);
            S=extra.assemb*KH*extra.assemb';
            % correct if direclet condtions is zero
            As   = S(~extra.disc.isdirdofs, ~extra.disc.isdirdofs);
            %% input proper fracture model
            problem.rock.perm(cells,:)=opt.max_perm*ones(nc,1);
            T = computeTrans(G, problem.rock);
            state_f = lincompTPFA(dt, state_f, G, T, pv, problem.fluid, problem.rock, 'MatrixOutput', true,...
                'wells', problem.W,'src', problem.src , 'bc', problem.bc_f);
            Af=state_f.A;
        end
        if(opt.do_plot && mod(count,opt.pstep)==0)
            figure(99)
            clf,%plotGrid(G),plotGrid(G,find(frac_cells),'FaceColor','r')
            val=sigm(:,2);
            val(frac_cells)=min(val(frac_cells));
            plotCellData(G,val)
            plotGrid(G,frac_cells,'FaceColor','k');colorbar
            if(sum(frac_cells)>0)
                if(strcmp(opt.frac_method,'mdem_pluss'))
                    eb=rlen>0 & e_broke(ihf);
                    plotFaces2D(G,hfn(~eb),'col','r','shift',hfshift(~eb))
                    plotFaces2D(G,hfn(eb),'col','b','shift',hfshift(eb))
                else
                    plotFaces2D(G,hfn(rlen<0),'col','r','shift',hfshift(rlen<0))
                    plotFaces2D(G,hfn(rlen>0),'col','b','shift',hfshift(rlen>0))
                end
            end
            pause(0.001)
        end
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
        not_converged=not(all(org_frac_cells == frac_cells) & all(org_e_broke ==e_broke));
        if(not_converged)
            %x=x_prev;
            count_it=count_it+1;
        end
    end
    states{count} = struct('t', t, 'pressure', p, 'uu', uu, 'bhp', x(ind_f(end)+1:end),'frac_cells',frac_cells,'count_it',count_it,'sigm',sigm,'stress',stress,'evec',evec);
    count = count+1;
    %{
    if(opt.do_plot && mod(count,opt.pstep)==0)
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
    %}
end
end




