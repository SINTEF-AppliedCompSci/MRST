% to do trap analysis we need  the following modules
%mrstModule add vertical-equil
%clear all
mrstModule add upscaling mimetic
%mrstModule('add','/home/hnil/heim/SVN/simmatlab/projects/statoil_upscaling')
%mrstModule('add','/home/hnil/heim/SVN/simmatlab/branches/halvor/mrst_extra/gridprocessing')
%addpath('/home/hnil/heim/SVN/simmatlab/branches/halvor/mrst_extra/gridprocessing')
theta_vec=[0 0.0162 0.03];
sat_mat=[];clear kr
gravity on;
clear Kkr
for j=1:numel(theta_vec)
    theta=theta_vec(j)

if(false)
    % we can load surfaces from the igems data set the possible names are
    % the directories  in data/surfaces 
    name='flatNP2';
    % there is 100 different surfaces we choose number
    i=2;
    % For testing it is usefull to use a coarse variant
    coarse=[4,4];
    %[Gt,rock2D,SVE,rock,G]=readIGEMSIRAP(name,i,coarse)
    [Gt,rock2D,rock,G]=readIGEMSIRAP(name,i,coarse)
else
    % simple test grid for finding traps
    L=0.3e3;L_p=L/2;H=50
    G=cartGrid([100 2 1],[L L H]);
    ff=@(x) sin(pi*x/L_p)*2+x*tan(theta);
    G.nodes.coords(:,3)=ff(G.nodes.coords(:,1))+G.nodes.coords(:,3);
    %100+G.nodes.coords(:,3)+G.nodes.coords(:,1)*0.0...
    %    -4*sin(pi*G.nodes.coords(:,1)/L_p);
    G=computeGeometry(G);
    Gt=topSurfaceGrid(G);
    % overide 2d surface calulalations
    Gt=computeGeometry(Gt);
    Gt.cells.z=ff(Gt.cells.centroids(:,1));
    Gt.faces.z=ff(Gt.faces.centroids(:,1));
    Gt.cells.H=ones(Gt.cells.num,1)*H;
    Gt.grav_pressure=@(g,omega) zeros(numel(g.cells.faces(:,1)),1);
end
clear G;
G=Gt;
clear Gt;
clear diff;
cellno=rldecode([1:G.cells.num]',diff(G.cells.facePos));
dhfz=-(G.cells.z(cellno)-G.faces.z(G.cells.faces(:,1)));
if(false)
    

N = G.faces.neighbors;
i = ~any(N==0, 2);
dzf=zeros(G.faces.num,1);
dzf(i)=G.cells.z(N(i,2)) - G.cells.z(N(i,1));
sgn=2*(G.faces.neighbors(G.cells.faces(:,1),1)==cellno)-1;
dzf_tmp=accumarray(G.cells.faces(:,1),dhfz.*sgn);
[dzf,dzf_tmp]
end
%return
%% plot traps
%%
%fluid = initSimpleVEFluidSForm('mu',[1 10],'rho',[700 1000],'sr',[0 0],'height',G.cells.H);
%fluid = initSimpleVEFluidSForm('mu',[6e-5 6.9e-4],'rho',[700 1020],'sr',[0 0],'height',G.cells.H);
fluid = initSimpleVEFluid_s('mu',[6e-5 6.9e-4],'rho',[760 1100],'sr',[0 0],'height',G.cells.H);
rock.poro=ones(G.cells.num,1);
rock.perm=ones(G.cells.num,1).*G.cells.H*darcy;
% define periodic boundary
dp_scale=1;
if(G.griddim>1)
   bcl{1}=pside([],G,'LEFT',0);dp{1}=0;
   bcr{1}=pside([],G,'RIGHT',0);
   bcl{2}=pside([],G,'BACK',0);dp{2}=0;
   bcr{2}=pside([],G,'FRONT',0);
else
   bcl{1}=struct('face',1,'value',dp_scale);dp{1}=dp_scale;
   bcr{1}=struct('face',G.faces.num,'value',dp_scale);
end
if(G.griddim>2)
   bcl{3}=pside([],G,'TOP',0);dp{3}=0;
   bcr{3}=pside([],G,'BOTTOM',0);
end
% find grid dimensions of the upscaling domain
L=nan(1,G.griddim);
for i=1:G.griddim
   L(i) = max(G.faces.centroids(:,i))-min(G.faces.centroids(:,i));
end

% make periodic grid
[Gp,bcp]=makePeriodicGridMulti3d(G,bcl,bcr,dp);
%[Gp,bcp]=makePeriodicGridMulti(G,bcl,bcr,dp);
% compute trans on the original grid
%used in transport
Trans = computeTrans(G,rock);
hfmap=Gp.cells.faces(:,end);
Transp=Trans(hfmap);
dhfzp=dhfz(hfmap);

%psolver_mim   =@(state,Gp,fluid_pure,bcp_new,rock) solveIncompFlow(state,Gp,computeMimeticIPGp(G,Gp,rock),fluid_pure,'bcp',bcp_new);
psolver_mim   =@(state,Gp,fluid_pure,bcp_new,rock) solveIncompFlow(state,Gp,computeMimeticIPGp(G,Gp,rock),fluid_pure,'bcp',bcp_new);
psolver=psolver_mim;
fluid_pure=initSingleFluid('mu',1,'rho',1);
perm = upscalePermeabilityPeriodic(Gp,bcp,1,psolver,fluid_pure,rock,L)
%return
%perm_mim = upscalePermeabilityPeriodic(Gp,bcp,dp_scale,psolver_mim,fluid_pure,rock,L);
%psolver_transport=@(state,Gp,fluid_nc,bcp) incompTPFA_extra(state,Gp,computeTransGp(G,Gp,rock),fluid_nc,'bcp',bcp,'trans_calc','upwind');
%psolver_transport=@(state,Gp,fluid_nc,bcp) incompTPFA(state,Gp,computeTransGp(G,Gp,rock),fluid,'bcp',bcp);
psolver_transport=@(state,Gp,fluid_nc,bcp) solveIncompFlow(state,Gp,computeMimeticIPGp(G,Gp,rock),fluid_nc,'bcp',bcp);
psolver_transport=@(state,Gp,fluid_nc,bcp) noflow(state);
% calculate porevolume
pv=poreVolume(G,rock);
% loop through all saturations
dS=1/100;
%sat_vec=[dS:dS:0.03];%1-dS];
sat_vec=linspace(0.000/50,4/50,50);%1-dS];
dp_vec=1;
sim_ok=false(numel(sat_vec),numel(dp_vec));
%{
Kkr=cell(2,1);
for kk=1:2;
 Kkr{kk}=nan(numel(sat_vec),numel(dp_vec),G.griddim^2);
end
kr=Kkr;
%}
%state=initResSol(Gp, 100*barsa, 0.1);%
   
% here capillary and visourse limit should be implemented   
%[pc_upsc, pc_max, pc_min, sat_min, sat_max] = upscalePc(G, fluid, pv);
%pc_values  = pc_upsc(sat_vec);
initcap=false;
dir=1;
bcp.value(bcp.tags==dir)=1;
state=initResSol(Gp, 100*barsa, 0.0);
state.s(:)=0.001;
state.extSat=[state.s,1-state.s];
fluid.s_min=0;
fluid.s_max=1;

for i=1:numel(sat_vec);      
  sat=sat_vec(i);
      if ~initcap
         %use previous state
         if true
            %sat=sat_vec(i);
            % define first guess for next stationary state
            if(false)
               eff_sat=sum(state.s.*pv)/sum(pv(:));
               if(eff_sat>0)
                  s_tmp=(state.s-eff_sat)+ones(size(state.s))*sat         ;
                  s_tmp=min(fluid.s_max-(1e-2),s_tmp);
                  s_tmp=max(fluid.s_min+(1e-2),s_tmp);
                  state.s=s_tmp;
                  assert(all(state.s(:)>=0))
                  assert(all(state.s(:)<=1))
               else
                  state.s=ones(size(state.s))*sat;
               end
            else
               state.s=ones(size(state.s))*sat;
            end
         end
      else
         % use capillary limit as guess of stationary state
         state.s=fluid.pcinv(repmat(pc_values(i),G.cells.num,1));
      end      
      
      %% search for stationary state
      mu=fluid.properties();
      V_i=perm(dir,dir)*dp_vec/(L(dir)*min(mu));
      DT_min=0.5*L(dir)/V_i
      
      %burde vært max?
      DT_min=min(DT_min,1*year);
      opt=struct( 'verbose',         false,...
         'nltol'   ,         1.0e-6, ...  % Non-linear residual tolerance
         'lstrials',         10    , ...  % Max no of line search trials
         'maxnewt' ,         20    , ...  % Max no. of NR iterations
         'tsref'   ,         1    , ...  % Time step refinement
         'max_it'  ,         100,...
         'max_newton',       200,...
         'DT_min',  DT_min,...
         'solve_pressure', true,...
         'dhfz',dhfzp);
     %[state,report] = simulateToSteadyStatePeriodic(state,Gp,bcp,psolver_extern,Transp,rock,fluid,fluid_nc,opt) 
      %[state,report] = simulateToSteadyStatePeriodic(state,Gp,bcp,psolver_transport,Transp,rock,fluid,[],opt);
      [state, report] = simulateToSteadyState(state, Gp, rock, fluid, DT_min,...
                         'psolver',psolver_transport,'bcp',bcp,'trans',Transp,...
                         'dhfz',dhfzp, 'verbose',false,...
                         'nltol',1.0e-6,'lstrials',10, 'maxnewt',20, 'tsref',1,'max_it',100,...
                         'max_newton',200,'solve_pressure', true);
      %%
      figure(1),clf
      plot(G.cells.centroids(:,1),G.cells.z,'b*',...
          G.cells.centroids(:,1),G.cells.z-state.s.*G.cells.H,'r*',...
          G.cells.centroids(:,1),G.cells.z-G.cells.H,'b*')      
      disp('hei')
      %%
      %pause
      %{
      figure(33)
      if(G.griddim>1)
         clf,plotCellData(G,state.s),colorbar
      else
         clf,
         zlevel=0.1*sin(3*2*pi*G.cells.centroids(:,1)/L(1));
         plot(G.cells.centroids(:,1),zlevel,'k')
         hold on
         plot(G.cells.centroids(:,1),zlevel-state.s,'b')
         hold off
      end
      drawnow;
      %}
      %pause
      %state=incompTPFA(state,Gp,Transp,fluid,'bcp',bcp);
      newton_steps(i,j)=report.newton_steps;
      steps(i,j)=report.steps;
      if(~report.stationary)
         sim_ok(i,j)=false;
      else
         sim_ok(i,j)=true;
      end
      % do single phase calculation to find phase fluxes
      kr_tmp=fluid.relperm(state.s,state);     
      Ts=cell(2,1);rocks=cell(2,1);
      for kk=1:2;
         % hack in case of non connected domains of fluid
         % should be fixed in pressure solver
         min_kr=min(kr_tmp((kr_tmp(:,kk)>0), kk));
         if isempty(min_kr)
            min_kr = 0;
         end
         
         %min_kr=min(kr_tmp(:,kk));
         aa=max(min_kr*1e-5,kr_tmp(:,kk));
         % multilply permeability with relperm
         rocks{kk}.perm=bsxfun(@times,rock.perm,aa);
         %only done to test to main
         Ts{kk} = computeTrans(G,rocks{kk});
         TT=Ts{kk}(Gp.cells.faces(:,end));
         % do single phase upscaling in all direction given the stationary
         % state
         if(~all(TT==0));
            Kkr_tmp=upscalePermeabilityPeriodic(Gp,bcp,dp_scale,psolver,fluid_pure,rocks{kk},L);
         else
            Kkr_tmp=zeros(G.griddim);
         end
         Kkr{kk}(i,j,:)=Kkr_tmp(:);
         aaa=Kkr_tmp/perm;
         kr{kk}(i,j,:)=aaa(:);
      end
      disp(['Saturation ',num2str(sat)])
      pv=poreVolume(G,rock);
      sat_mat(i,j)=sum(pv.*state.s(:,1))/sum(pv);
end
end
%%
figure(2),clf,hold on
plot(sat_mat*50,kr{1}(:,:,1),'*-')
krCO2=kr{1};
krW=kr{2};
save('data/upscaled_relperm_theta','sat_mat','krCO2','krW')


