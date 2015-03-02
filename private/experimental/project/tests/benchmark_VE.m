
Lx=1e4;Ly=1e4;
K=(100e-3)*darcy();
phi= 0.1;
dt=0.1*year;
rate = 1e6;
thread_vec=[0,2.^[0:3]];
 wtime_vec=[]
 nv_vec=[];
for nn=1:4:20;
   g=cartGrid([10 10 1]*nn,[nn*1000 nn*1000 10]);
   g=twister(g);
   g.nodes.coords(:,3)=g.nodes.coords(:,3)...
      +sin(g.nodes.coords(:,1).*2*pi/Lx).*sin(g.nodes.coords(:,2).*2*pi/Ly);
   g=computeGeometry(g);
   g_top = topSurfaceGrid(g);
   fluid = initVEFluid(g_top, 'mu', [0.1 0.4]*centi*poise, 'rho', [600 1000], 'sr', 0, 'sw', 0);
   rock.perm = K*ones(g.cells.num,1);
   rock.poro = phi*ones(g.cells.num,1);
   solinit = initResSol(g_top, 0);
   solinit.h=0*solinit.s;
   solinit.max_h=solinit.h;   
   st=tic;
   bc_faces=find(sum(g_top.faces.neighbors>0,2)==1);
   bc = addBC([], bc_faces, 'pressure', g_top.faces.z(bc_faces)*1000*norm(gravity));
   bc.h=0;
   cind=floor(g_top.cartDims/2)+1;
   ind_well=cart2active(g_top,sub2ind(g_top.cartDims,cind(1),cind(2)));
   W = addWell([], g_top, rock, ind_well,'Type', 'rate','Val',rate/year,'Radius',0.1,'Dir','z');
%% use MRST to add wells we need to change defintion of well
   disp(['Testing grid of size ',num2str(g_top.cells.num/1e3) ' kc'])
for i=1:numel(W)
   W(i).compi=nan;
   W(i).h = g_top.cells.H(W(i).cells);
end
   S_2d = computeMimeticIPVE(g_top, rock,'Innerproduct','ip_simple');
   solinit = solveIncompFlowVE(solinit, g_top, S_2d, rock, fluid, 'bc', bc,'wells',W);
   wtime=zeros(numel(thread_vec),1);
   st=tic;
   sol_ve = explicitTransportVE(solinit,g_top, dt, rock, fluid,'time_stepping','simple','Verbose',false,'bc',bc,'wells',W);
   wtime(1)=toc(st);
   for i=1:numel(thread_vec);
      nt=thread_vec(i);
      VETransportCPU();
      sol_va=solinit;
      st=tic;
      [sol_va.h,sol_va.max_h] = VETransportCPU(solinit, g_top, dt, rock, fluid, 'computedt',...
        true,'verbose',false,'gravity', norm(gravity), 'flag', nt,'wells',W); 
      wtime(1+i)=toc(st);
   end
   wtime_vec=[wtime_vec,wtime];
   nv_vec=[nv_vec,g_top.cells.num];
end
%%
leg={'VA'}
for nt=thread_vec;
   leg{end+1}=['VE', num2str(nt)];
end
figure(1),plot(nv_vec/1e3,wtime_vec,'*-')
legend(leg{:})
figure(2),clf,
subplot(2,1,1)
plotCellData(g_top,sol_va.h);caxis([0,10])
subplot(2,1,2)
plotCellData(g_top,sol_ve.h);caxis([0,10])
%%
%
figure(3),plot(nv_vec/1e3,1./bsxfun(@rdivide,wtime_vec(3:end,:),wtime_vec(3,:)),'*-')
legend(leg{3:end})
%figure(3),plot(nv_vec/1e3,1./bsxfun(@rdivide,wtime_vec(2:end,:),wtime_vec(2,:)))
