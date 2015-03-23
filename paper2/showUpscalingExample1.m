%% plot saved results from stored from upscalingExample1
upAquifer = makeAquiferModel('A',0);
fAquifer  = makeAquiferModel('A',2);

z  = fAquifer.Gt.cells.z;
zt = max(z)*ones(size(z));
for i=2:numel(zt)-1
   zt(i)=max(z(i:end));
end
zt(end) = max(zt(end-1),z(end));
ht  = zt - z;
ff = exp(-linspace(-25,25,501).^2); ff=ff'/sum(ff);
hts = filter2(ff, ht);

legendtext = {'Fine scale', 'Accretion layer', ...
              'Analytic: square', 'Analytic: sinus'};
linetype = {'k-', 'b-','r-','g-'};

surf_topos = {'smooth','inf_rough','square','sinus'};
fluid = makeFluidModel(upAquifer, 'residual', true, ...
         'dissolution', false, 'fluidType', 'sharp interface', ...
         'top_trap', hts, 'surf_topo',  surf_topos{1});
res=load('data/upscalingExample1Data.mat');
results=res.results;
for i=1:numel(surf_topos)   
   pstep=numel(res.control.step.val)-70;% choose time step from stored data
   fprintf('Plot time step %4i year for %s \n', floor(sum(res.control.step.val(1:pstep))/year), legendtext{i});
   state = results{i}.states{pstep};% 
   
   sG = free_sg(state.s(:,2),state.smax(:,2), ...
      struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));% calulate the free CO2
   hold on
   plot(res.xc, filter2(results{i}.ff,sG.*upAquifer.Gt.columns.dz), linetype{i}, 'LineWidth', 2);     
end
axis tight
set(gca,'YDir','reverse','FontSize',12);
h=legend(legendtext{:}, 4); set(h,'FontSize',14);
set(gcf,'Position', [680 580 800 420]);