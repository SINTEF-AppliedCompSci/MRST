do_print=true;
if(do_print)
 mkdir('figs')
end
gravity on;
gravity on;
depth=1300;
res=load(['data/disolutionExample1Data_',num2str(depth),'.mat']);
legendtext = {'No dissolution (A=0)', 'Dissolution (A=0)', ...
    'No dissolution (A=2)', 'Dissolution (A=2)'};
linetype = {'b-', 'b--', 'r-', 'r--'};
upAquifer = makeAquiferModel('A',0,'D',depth);
fAquifer  = makeAquiferModel('A',2,'D',depth);



z  = fAquifer.Gt.cells.z;
zt = max(z)*ones(size(z));
for i=2:numel(zt)-1
   zt(i)=max(z(i:end));
end
zt(end) = max(zt(end-1),z(end));
ht  = zt - z;
ff = exp(-linspace(-25,25,501).^2); ff=ff'/sum(ff);
hts = filter2(ff, ht);


% left and right panel of figure 15
%%


residual_c=true;% plot residual case
n_c=1;% plot top surface case
dissolution_c=true%plot disolution case

% find number in res
k=1;
for residual= [false,true] %residual saturation or not   
    for n=1:2, % flat or non flat topsurface
        %% Make fluid model
        for dissolution=[false,true]
            if(n_c==n && dissolution_c==dissolution  && residual_c==residual)
                %break;
                kk=k;
            end            
            k=k+1;
        end
    end
end
if(n_c==1)
  aquifer=upAquifer;  
else
  aquifer=fAquifer;  
end

if(n_c==1)
   h_trap=0; 
else
    h_trap=ht;
end
xc=res.xc;
zt=aquifer.Gt.cells.z*0;
zb=aquifer.Gt.cells.z*0+aquifer.Gt.cells.H;
Gt=aquifer.Gt;
for nn=[40,50,70]
    figure(nn),clf
    fprintf('Year from start injection %i year\n', floor(sum(res.control.step.val(1:nn-1))/year))
    fprintf('Year after stop of injection %i year\n', floor(sum(res.control.step.val(1:nn-1))/year)-50)
    state = res.results{kk}.states{nn};
    fluid = makeFluidModel(aquifer, 'residual', residual, ...
        'dissolution', dissolution, 'fluidType', 'sharp interface','only_pvt',false);%,'co2_type','coolprops');
    if(~dissolution)
        sG_free = free_sg(state.s(:,2),state.smax(:,2), ...
            struct('res_gas',fluid.res_gas, 'res_oil', fluid.res_oil));
        sGmax=state.smax(:,2);
    else
        sG_free = free_sg(state.s(:,2),state.sGmax, ...
            struct('res_gas',fluid.res_gas, 'res_oil', fluid.res_oil));
        sGmax=state.sGmax;
    end
    
    %%
    p=state.pressure;sG=state.s(:,2);
    rs=state.rs;   
    diff=max(zt)-min(zt); 
    
    
    h=(sG_free.*Gt.cells.H)./(1-fluid.res_oil);
    
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
      [zt + h_trap; zt(end:-1:1)], myCOColor(1))
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_trap(end:-1:1)], myCOColor(2))
 
    h_max=(sGmax.*Gt.cells.H)./(1-fluid.res_oil);
    if(dissolution)
        mm=minRs(p,sG,sGmax,fluid,Gt);%.*Gt.cells.H;
        h_res=Gt.cells.H.*((1-sG).*state.rs-mm)/fluid.dis_max;
    else
        h_res=0*Gt.cells.H;
    end
    
        
    assert(all(h_res>-1e-4))
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_max(end:-1:1)],myCOColor(3))
    patch(xc([1:end end:-1:1]), ...
      [zt + h_max; zt(end:-1:1)+h_max(end:-1:1)+h_res(end:-1:1)], myCOColor(4))
    patch(xc([1:end end:-1:1]), ...
      [zt + h_max+h_res; zb(end:-1:1)], myCOColor(5))
  
    line(xc([1:end]),zt + h_max,'LineWidth',2,'Color','k');%#ok
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    set(gca,'YDir','reverse'), axis tight
    %%
    axis([0 25 0 50])
    if(do_print)
        %print('-depsc2',[fname,'_fig_',num2str(nn),'panel.eps'])
        print('-depsc2',['figs/fig14_',num2str(nn),'_panel.eps'])
    end
    drawnow;   
    pause(0.1)
end

    
    
    

