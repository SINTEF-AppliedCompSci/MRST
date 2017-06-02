%% 1D Aquifer: numerical upscaling of relative permeabilities
% Compute numerically upscaled relative permeabilities for a set of
% representative dip angles.

mrstModule add co2lab upscaling incomp mimetic;
gravity on;

% Preallocate memory
theta_vec=[0 0.0162 0.03];
N = 51;
n = numel(theta_vec);
sat_mat=nan(N,n);
t = cell(1,2); t{1}=zeros(N,n,4); t{2}=zeros(N,n,4);
[kr,Kkr] = deal(t); clear t;

hf = figure;

for j=1:numel(theta_vec)
   theta=theta_vec(j);

   %% Create grid and fluid models
   [L,H] = deal(.3e3, 50);
   G = cartGrid([100 2 1],[L L H]);
   ff = @(x) sin(2*pi*x/L)*2+x*tan(theta);
   G.nodes.coords(:,3) = ff(G.nodes.coords(:,1)) + G.nodes.coords(:,3);
   G = topSurfaceGrid(computeGeometry(G));
   G = computeGeometry(G);

   G.grav_pressure = @(g,omega) zeros(numel(g.cells.faces(:,1)),1);
   cellno = rldecode((1:G.cells.num)',diff(G.cells.facePos));
   dhfz   = -(G.cells.z(cellno)-G.faces.z(G.cells.faces(:,1)));

   fluid = initSimpleVEFluid_s( 'mu',[6e-5 8e-4], 'rho',[760 1100], ...
      'sr',[0 0],'height', G.cells.H);
   rock.poro = ones(G.cells.num,1)*.2;
   rock.perm = ones(G.cells.num,1)*100*milli*darcy;

   %% Setup periodic grid
   [bcl,bcr,dp]=deal(cell(2,1));
   bcl{1} = pside([],G,'LEFT', 0);  dp{1}=0;
   bcr{1} = pside([],G,'RIGHT',0);
   bcl{2} = pside([],G,'BACK', 0);  dp{2}=0;
   bcr{2} = pside([],G,'FRONT',0);

   L = nan(1,G.griddim);
   for i=1:G.griddim
      L(i) = max(G.faces.centroids(:,i))-min(G.faces.centroids(:,i));
   end
   [Gp,bcp] = makePeriodicGridMulti3d(G,bcl,bcr,dp);
   bcp.value(bcp.tags==1) = 1;

   %% Setup solvers
   pv    = poreVolume(G,rock);
   hT    = computeTrans(G,rock);
   hfmap = Gp.cells.faces(:,end);
   hTp   = hT(hfmap);
   dhfzp = dhfz(hfmap);
   psolver = @(state, Gp, fluid_pure, bcp_new, rock) ...
      incompMimetic(state, Gp, computeMimeticIPGp(G,Gp,rock), ...
                    fluid_pure, 'bcp', bcp_new);
   fluid_pure = initSingleFluid('mu',1,'rho',1);
   perm = upscalePermeabilityPeriodic(Gp,bcp,1,psolver,fluid_pure,rock,L);
   psolver_transport= @(state,Gp,fluid_nc,bcp) noflow(state);

   state        = initResSol(Gp, 100*barsa, 0.0);
   state.s(:)   = 0.001;
   state.extSat = [state.s,1-state.s];
   fluid.s_min  = 0;
   fluid.s_max  = 1;

   %% Main loop
   fprintf(1,'\nUpscaling relative permeability for angle=%f\n', theta);
   fprintf(1,'Saturation values %f to %f\n', 0, 4/H);
   sat_vec = linspace(0.0,4/H,N);
   sim_ok  = false(numel(sat_vec),1);
   for i=1:numel(sat_vec);
      sat = sat_vec(i);
      state.s = ones(size(state.s))*sat;

      %% search for stationary state
      mu     = fluid.properties();
      V_i    = perm(1,1)*1.0/(L(1)*min(mu));
      DT_min = min(0.5*L(1)/V_i, year);

      [state, report] = simulateToSteadyState(state, Gp, rock, fluid, DT_min,...
         'psolver',psolver_transport,'bcp',bcp,'trans',hTp,...
         'dhfz',dhfzp, 'verbose',false,...
         'nltol',1.0e-6,'lstrials',10, 'maxnewt',20, 'tsref',1,'max_it',100,...
         'max_newton',200,'solve_pressure', true);

      figure(hf); clf;
      plot(G.cells.centroids(:,1), G.cells.z,'b*',...
         G.cells.centroids(:,1), G.cells.z-state.s.*G.cells.H,'r*',...
         G.cells.centroids(:,1), G.cells.z-G.cells.H,'b*')
      if(~report.stationary)
         sim_ok(i,j)=false;
      else
         sim_ok(i,j)=true;
      end
      drawnow;
      %% Single phase calculation to find phase fluxes
      kr_tmp = fluid.relperm(state.s,state);
      for kk=1:2;
         % hack in case of non connected domains of fluid
         % should be fixed in pressure solver
         min_kr = min(kr_tmp(kr_tmp(:,kk)>0, kk));
         if isempty(min_kr)
            min_kr = 0;
         end

         aa = max(min_kr*1e-5, kr_tmp(:,kk));
         srock.perm = bsxfun(@times, rock.perm, aa);
         hT = computeTrans(G, srock);
         hT = hT(Gp.cells.faces(:,end));

         % Single-phase upscaling in all direction given the stationary
         % state
         if(~all(hT==0));
            Kkr_tmp = upscalePermeabilityPeriodic(Gp, bcp, 1,...
               psolver, fluid_pure, srock, L);
            Kkr{kk}(i,j,:) = Kkr_tmp(:);
            kr_t           = Kkr_tmp/perm;
            kr{kk}(i,j,:)  = kr_t(:);
         else
            tmp = zeros(G.griddim);
            Kkr{kk}(i,j,:) = tmp(:);
            kr{kk}(i,j,:)  = tmp(:);
         end
      end
      if ~mod(i,10)
         fprintf(1,' %.3f\n',sat);
      else
         fprintf(1,' %.3f',sat);
      end
      sat_mat(i,j)=sum(pv.*state.s(:,1))/sum(pv);
      %%
   end
end

%% Plot final results
figure,clf, hold on
plot(sat_mat*H,kr{1}(:,:,1),'*-')
krCO2=kr{1};
krW=kr{2};
filename = fullfile(mrstPath('co2lab'), 'examples', 'papers', 'COMG-1', 'data', 'upscaled_relperm_theta');

ensure_path_exists(filename);
save(filename,'sat_mat','krCO2','krW')
