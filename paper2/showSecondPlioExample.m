%% Inspect the Pliocenesand Formation
% This formation has a very low percentage (0.02%) of trapping compared to
% the overall volume of the whole model
try
    require coarsegrid deckformat mex ad-fi ad-props
catch %#ok<CTCH>
    mrstModule add coarsegrid deckformat mex ad-fi ad-props
end
do_print=true;
grdecl = getAtlasGrid('Pliocenesand');
G      = processGRDECL(grdecl{1});
G      = computeGeometry(G(1));
Gt     = topSurfaceGrid(G);
%ta     = trapAnalysis(Gt, false);



%%
% Next, we will perform a VE simulation. This could, of course, have been
% launched from inside the interactive viewer, but to make the example as
% reproducible as possible, we launch it manually from the outside.

%%
% Rerun for a longer time
petrodata.avgperm = 1.2*darcy;
petrodata.avgporo = 0.25;

%%
% Adding depth to the plio example
% It is to shallow for a real storage sight
depth=1200;


% cut grid to avoid calculation on not relevant domain
wpos=Gt.parent.cells.centroids(5280,1:2);
wpos(:,1)=4.85e5;
G=Gt.parent;
rm_cells=abs(Gt.cells.centroids(:,2)-wpos(:,2))>2.5e4;
G=removeCells(G,rm_cells);
G.nodes.coords(:,3)=G.nodes.coords(:,3)+depth;
G=computeGeometry(G);
Gt=topSurfaceGrid(G);
ta     = trapAnalysis(Gt, false);
%clear G;

% set the permeability and porosity
rock.poro = repmat(petrodata.avgporo, Gt.cells.num, 1);
rock.perm = repmat(petrodata.avgperm, Gt.cells.num, 1);
rock2D    = averageRock(rock, Gt);
%% Find pressure boundary
% Setting boundary conditions is unfortunately a manual process and may
% require some fiddling with indices, as shown in the code below. Here, we
% need to find all outer vertical faces
i = any(Gt.faces.neighbors==0, 2);  % find all outer faces
I = i(Gt.cells.faces(:,1));         % vector of all faces of all cells, true if outer
j = false(6,1);                     % mask, cells can at most have 6 faces,
%j(1:4)=true;
j(1)=true;
%   extract east, west, north, south
J = j(Gt.cells.faces(:,2));         % vector of faces per cell, true if E,W,N,S
bcIxVE = Gt.cells.faces(I & J, 1);

%% Set time and fluid parameters

% Fluid data are taken from paper SPE 134891
gravity on


%%
%bc = addBC([], bcIxVE, 'pressure', ...
%    Gt.faces.z(bcIxVE)*rhow*norm(gravity));
%bc.sat = zeros(size(bc.face));

aquifer.G  = G;
aquifer.Gt = Gt;
aquifer.rock   = rock;
aquifer.rock2D = rock2D;
%aquifer.W  = W;
%aquifer.bc = bc;

residual=true;
res={};
for k=1:3
    switch k
        case 1
            dissolution=false;
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface');
            %system = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
            %    fluid, 'simComponents', s, 'VE', true);
        case 2
            dissolution=true;
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface');
            fluid=rmfield(fluid,'dis_rate');
            %system = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
            %    fluid,'simComponents',s,'VE',true);
        case 3
            % deside which fluid to use
            dissolution=true;
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface');
            %system = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
            %    fluid,'simComponents',s,'VE',true);
        otherwise
            error('No such case')
    end
    
    z  = aquifer.Gt.cells.z;
    %plotting
    % Plot contours
    %{
    if(dissolution)
        systemOG = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
            fluid,'simComponents',s,'VE',true);
    else
        systemOG = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
            fluid, 'simComponents', s, 'VE', true);
    end
    systemOG.nonlinear.linesearch    = false;
    systemOG.nonlinear.maxIterations = 10;
    systemOG.nonlinear.tol           = 1e-6;
    t2=tic;
    [wellSols, states] = runMrstADI(x0, Gt, system, control, ...
        'force_step', false, 'dt_min', 0.5*year, 'report_all', false, 'bc', bc);
    t2=toc(t2);
    %}
    res{k}=load(['data/secondPlioExample_',num2str(depth),'_',num2str(k),'.mat']);%#ok   
    fluids{k}=fluid;%#ok
end


%% Pliocenesand: effect of dissolution
Gt = aquifer.Gt;

leg_traps = {'Dissolved','Residual (traps)', 'Residual', ...
   'Residual (plume)', 'Movable (traps)', ...
   'Movable (plume)', 'Leaked'};
leg_methods = {'Compressible', ...
   'Instant dissolution','Dissolution rate', ''};

% Plot contours
methods   = [1,3,2];
mysteps   = [69 91];% time 940 1490
%ptimes=[950 1500]
%ptimes=[955 1455]
%ptimes=[955 1455]
ptimes=[710 1510];

%ptimes=[100 1500]
%mysteps   = [10 30];% time 940 1490
leg_methods=leg_methods(methods);
m = numel(methods);
col = .5*(lines(m) + ones(m,3));
x = Gt.cells.centroids(:,1);
y = Gt.cells.centroids(:,2);
W=res{1}.control.W{1};
time=[0;cumsum(res{1}.control.step.val)/year];
for k=1:numel(ptimes)
   [d,i]=min((time-ptimes(k)).^2);
   mysteps(k)=i;
end


ptimes=time(mysteps);
for step=mysteps
   figure, hold on;
   has_rs=[];
   for n=1:3     
       time_tmp=[0;cumsum(res{methods(n)}.control.step.val)/year];
       assert(time_tmp(step)==time(step));
      % calculate properties
      state=res{methods(n)}.states{step};
      p= state.pressure;
      sG= state.s(:,2);
      if(n>1)
          sGmax=state.sGmax;
      else
          sGmax=state.smax(:,2);
      end
      sG = free_sg(sG, sGmax, ...
                    struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
      fluid=fluids{methods(m)};
      drho=fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p); 
      h=(fluid.pcWG(sG,p,'sGmax',sGmax))./(drho*norm(gravity()));
      h_max=(fluid.pcWG(sGmax,p,'sGmax',sGmax))./(drho*norm(gravity()));
      if(n==3)
          mm=minRs(p,sG,sGmax,fluid,Gt);%.*Gt.cells.H;
          h_res_diff=Gt.cells.H.*((1-sG).*state.rs-mm)/fluid.dis_max;
      else
            h_res_diff=0*Gt.cells.H;
      end
      %h_res_diff=((state.rs.*state.s(:,1)-min_rs)./fluid.dis_max).*Gt.cells.H;
      assert(all(h_res_diff>=-1e-5));
      h_res=h_max+h_res_diff;  
      % calulate masses for later
%{
      masses = phaseMassesVEADI(Gt, state, rock2D, fluid);
      co2mass= masses(1)+masses(3);
      %% h and h_max;
      state.h=h;
      state.h_max = h_max;
      % calculated distributions only valid for sharp interface.
      masses = massTrappingDistributionVEADI(Gt, state, rock2D, fluidADI, ts);
%}      
      
      FF_cc = TriScatteredInterp(x, y, h_max);%#ok
      [xx,yy] = meshgrid(linspace(min(x),max(x),200),...
                         linspace(min(y),max(y),200));
      cc=FF_cc(xx,yy);
      cc(cc<.1)=NaN;
      [c,a] = contourf(xx, yy, cc, 1e-1);
      set(get(a,'Children'),'FaceColor',col(n,:),'EdgeColor',col(n,:));
      set(a,'EdgeColor',col(n,:));
   end
   FF_zz = TriScatteredInterp(x,y,Gt.cells.z);%#ok
   zvec = FF_zz(xx,yy);
   contour(xx,yy,zvec,40,'k');
      
   axis equal tight off;   
   plot(Gt.cells.centroids(W.cells,1), ...
      Gt.cells.centroids(W.cells,2),'ok', ...
      'MarkerSize',10,'MarkerFaceColor','black')
   [legh,objh,outh,outm] = legend(leg_methods,4);
   set(legh,'Fontsize',12);
   for n=m+1:2*m
      set(get(get(objh(n),'Children'),'Children'),'LineWidth',4);
   end
   title(['after ' num2str(floor(time(step))),' years'],'FontSize',12)
   if(do_print)
      print('-depsc2',['figs/plio-large-1200-',num2str(floor(time(step))),'.eps'])
   end
end

%% calculating alle masses
%   has_rs=[];
%res={};
for m=1:3
    nsteps=numel(res{m}.states);
    res{m}.tot_masses=nan(nsteps,3);%#ok
    res{m}.masses=nan(nsteps,8);%#ok
    control=res{m}.control;
    for step=1:numel(res{m}.states)
        % calculate properties
        state=res{m}.states{step};
        p= state.pressure;
        sG= state.s(:,2);
        if(m==3)
            sGmax=state.sGmax;
        else
            sGmax=state.smax(:,2);
        end
        sG_free = free_sg(sG, sGmax, ...
            struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
        if(m==2)
            sGmax=-((1-fluid.res_water-fluid.res_gas)*sG_free -(1 - fluid.res_water) * sG)./fluid.res_gas;
        end
        assert(all(sGmax>=0));
        
        fluid=fluids{m};
        drho=fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p);
        h=(fluid.pcWG(sG,p,'sGmax',sGmax))./(drho*norm(gravity()));
        h_max=(fluid.pcWG(sGmax,p,'sGmax',sGmax))./(drho*norm(gravity()));
        if(n==3)
            mm=minRs(p,sG,sGmax,fluid,Gt);%.*Gt.cells.H;
            h_res_diff=Gt.cells.H.*((1-sG).*state.rs-mm)/fluid.dis_max;
        else
            h_res_diff=0*Gt.cells.H;
        end
        h_max_tmp=(sGmax/(1-fluid.res_water)).*Gt.cells.H;
        h_tmp=(sG_free/(1-fluid.res_water)).*Gt.cells.H;
        
        ss=(h*(1-fluid.res_water)+(h_max-h)*fluid.res_gas)./Gt.cells.H;
        ss_max=((1-fluid.res_water)*h_max)./Gt.cells.H;
        assert(all(abs(state.s(:,2)-ss)<1e-3))
        assert(all(abs(h_max-h_max_tmp)<1e-3));% valid for sharp interface when mass calculation is valid
        assert(all(abs(h-h_tmp)<1e-3));% valid for sharp interface when mass calculation is valid
        
        assert(all(h_res_diff>=-1e-5));
        h_res=h_max+h_res_diff;
        % calulate masses for later
        if(step==1)
          totMass=0;
        else
            if(control.step.control(step-1)==1)
                dT=control.step.val(step-1);
                totMass = totMass + fluid.rhoGS*W(1).val*dT;
            else                                
                % is no well
            end
        end
        
        masses = phaseMassesVEADI(Gt, state, rock2D, fluid);
        res{m}.tot_masses(step,:)=masses;%#ok
        co2mass= masses(1)+masses(3);
        % h and h_max;
        state.h=h;
        state.h_max = h_max;
        % calculated distributions only valid for sharp interface.
        dh=Gt.cells.z*0;%no subscale trapping
        masses = massTrappingDistributionVEADI(Gt, state, rock2D, fluid, ta,dh);
        % store all masses as a matrix with rows repersenting time
        res{m}.masses(step,:)=[masses,totMass];%#ok
    end
end



%% Draw inventory plot

figure(); set(gcf,'Position',[100 100 1400 600]);
xw = 0.8/numel(methods);
xa = 0.2/(numel(methods)+1);
for i=1:numel(methods)
   ind=[1:5,7,8]; 
   masses=res{methods(i)}.masses(:,ind);   
   masses(:,end) = masses(:,end)-sum(masses(:,1:end-1),2);      
   axes('position',[(i-1)*xw+i*xa, .17, xw, .78]);%#ok
   set(gca,'FontSize',14);
   hp = area(time, masses/1e9);
   col = getInventoryColors(1:7);
   for j=1:numel(hp), set(hp(j),'FaceColor',col(j,:)); end
   %xlabel('Years since simulation start');
   %ylabel('Mass (Mtonn)');
   axis tight
   title(leg_methods{i});
   %axis([0 1500 0 260])
   
end
axes('position',[0.05 0.01 0.9 0.03],'visible','off');
hL=legend(hp,leg_traps);
set(hL,'orientation', 'horizontal','position',[0.05 0.02 0.9 0.03]);
if(do_print)
    set(gcf,'PaperPositionMode','auto')
      print('-depsc2',['figs/plio-large-1200-inventory.eps']);%#ok
end


