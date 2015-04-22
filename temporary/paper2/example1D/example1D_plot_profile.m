%%
smooth=false;
depth=1300;
res_fluid=true;
mydir='test_data_dis_all/';fluid.dis_max=0.03; 
%mydir='test_data_dis_linear/';fluid.dis_max=0.02; 
mydir='data/';
use_dis=true;
if(smooth)
  if(res_fluid)
            fname=['smooth_res_fluid','_','_depth_',num2str(depth)]; 
        else
            fname=['smooth_nores_fluid','_','_depth_',num2str(depth)];
        end
    else
        if(res_fluid)
            fname=['res_fluid','_','_depth_',num2str(depth)];
        else
            fname=['nores_fluid','_','_depth_',num2str(depth)];
        end
end
if(use_dis)
    fname=[fname,'_dis'];
end   
res=load([mydir,fname,'.mat']);
   
   %%
%% Plot results
%figure
states=res.results{1}.states;
%p0=x0.pressure;
xc = res.xc;
zt = res.z*0;
zb = zt + 50;
mass=nan(numel(states),2);
if(res_fluid)
    sr= 0.21;sw= 0.11;%kwm= [0.75 0.54];
else
    sr=0;sw=0;kwm=[1 1];
    %sr= 0.21
end
fluid.res_gas=sr;
fluid.res_water=sw;
opt.res_gas=sr;
opt.res_water=sw;
example1D_makeparam
%if(use_dis)
%   fluid.dis_max=0.03; 
%end
%fluid=fluid{1};

%for nn=1:numel(states)
for nn=[40,50,70]    
    nnn=nn
    %nnn=20;
    figure(1),clf
    state=states{nn};
    p=state.pressure;sG=state.s(:,2);sGmax=state.sGmax;
    rs=state.rs;
    bW=fluid.bW(p,rs);bG=fluid.bG(p);
    pv=fluid.pvMultR(p).*Gt.cells.volumes.*Gt.cells.H.*rock.poro;
    mass(nn,:)=[sum(pv.*(bG.*sG + rs.*bW.*(1-sG))),sum(pv.*(rs.*bW.*(1-sG)))];%,sum(pv.*(bG.*))) ; 
    if(isfield(fluid,'dis_max'))
        smax=state.sGmax;
        sGmax=state.sGmax;
        sG=state.s(:,2);
        %rsH=Gt.cells.H.*state.s(:,1).*state.rs/fluid.dis_max;
        rsH=Gt.cells.H.*(state.s(:,1).*state.rs-...
            (sGmax-sG)*fluid.dis_max*(1-fluid.res_gas)/(1-fluid.res_water)-...
            (sG)*fluid.dis_max*(fluid.res_water)/(1-fluid.res_water) )/fluid.dis_max;
        assert(all(rsH>-1e-4))
        rsH=max(rsH,0)
     else        
        smax=state.smax(:,2);
        rsH=state.s(:,1)*0;
     end
    diff=max(zt)-min(zt); 
    
    
    % plot as VE
    subplot(3,1,2),cla   
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 300]);
    subplot(3,1,3),cla   
    plot(xc,(state.pressure-p0)/barsa); %set(gca,'YLim',[0 200]);
    subplot(3,1,1),cla,hold on
    if(isfield(fluid,'dis_max'))
    plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H,xc,smax.*Gt.cells.H+rsH); %set(gca,'YLim',[0 20]);
    else
       plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H);
       hh=state.s(:,2).*Gt.cells.H;
       plot(xc,filter2(ff./sum(ff),hh))
    end
    %%
    %zt=zt*0;
    figure(nnn)
     clf
     sG=free_sg(state.s(:,2),smax,opt);
    %
     h=(sG.*Gt.cells.H)./(1-fluid.res_water);
    %c(7)=patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    %c(6)=patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    clear c;
    c(4)=patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)], myCOColor(2))
  %
  %%{
    %sGmax=state.sGmax;
   
    h_max=(sGmax.*Gt.cells.H)./(1-fluid.res_water);
    h_res=(-((sG.*fluid.res_water/(1-fluid.res_water)+(sGmax-sG)*(1-fluid.res_gas)/(1-fluid.res_water))*fluid.dis_max-state.s(:,1).*state.rs)).*Gt.cells.H./fluid.dis_max;
    assert(all(h_res>-1e-4))
    c(1)=patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_max(end:-1:1)],myCOColor(3))
    c(2)=patch(xc([1:end end:-1:1]), ...
      [zt + h_max; zt(end:-1:1)+h_max(end:-1:1)+h_res(end:-1:1)], myCOColor(4))
    c(3)=patch(xc([1:end end:-1:1]), ...
      [zt + h_max+h_res; zb(end:-1:1)], myCOColor(5))
  %}  
    %line(xc([1:end]),zt + h_max,'LineWidth',2,'Color','k')
    %c(4)=patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20], .7*[1 1 1]);
    for kk=1:numel(c)
        set(c(kk),'EdgeColor','none')
    end
    set(gca,'YDir','reverse'), axis tight
    %%
    axis([0 30 0 50])
    print('-depsc2',['figs/',fname,'_fig_',num2str(nn),'panel.eps'])
    drawnow;
    
    pause(0.1)
end