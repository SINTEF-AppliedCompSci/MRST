%%
addpath(fullfile(pwd,'..'))
smooth=true;
for smooth=[true,false]
%for smooth=[false]
depth=2300;
res_fluid=true;
mydir='test_data_dis_all/';fluid.dis_max=0.03;
mydir='data_all_dis/';fluid.dis_max=0.03;
mydir='data_all_results_tol_1e_4/';fluid.dis_max=0.03;
mydir='data_all_changed_1e_6/';dis_max=0.03
mydir='data/';
%mydir='test_data_dis_linear/';fluid.dis_max=0.02; 
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
%for smooth=[true]%false
%surf_topos={'smooth','square','inf_rough','sinus'}
for mf=1:3
%mf=1;
if(~smooth)
    surf_topo='inf rough';
else
    surf_topo='smooth';
end
ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
%fluid_case='linear cap.';
fluid_case=ff_names{mf};
if(res_fluid)
    %
end
opt=res.opt;
example1D_makeparam
states=res.results{mf}.states;
%p0=x0.pressure;
xc = res.xc;
zt = res.z;
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

%if(use_dis)
%   fluid.dis_max=0.03; 
%end
%fluid=fluid{1};

%for nn=1:numel(states)
for nn=[50,80]    
    nnn=nn
    %nnn=20;
    %figure(1),clf
    state=states{nn};
    %p=state.pressure;sG=state.s(:,2);sGmax=state.sGmax;
    %rs=state.rs;
    %bW=fluid.bW(p,rs);bG=fluid.bG(p);
    %pv=fluid.pvMultR(p).*Gt.cells.volumes.*Gt.cells.H.*rock.poro;
    %mass(nn,:)=[sum(pv.*(bG.*sG + rs.*bW.*(1-sG))),sum(pv.*(rs.*bW.*(1-sG)))];%,sum(pv.*(bG.*))) ; 
    if(isfield(fluid,'dis_max'))
        smax=state.sGmax;
        sGmax=state.sGmax;
        sG=state.s(:,2);
        press=state.pressure;
        min_rs=   minRs(press,sG,sGmax,fluid, Gt);
        assert(all((state.rs-min_rs./state.s(:,1)))>-1e-5);
        %rsH=Gt.cells.H.*state.s(:,1).*state.rs/fluid.dis_max;
        %{
        rsH=Gt.cells.H.*(state.s(:,1).*state.rs-...
            (sGmax-sG)*fluid.dis_max*(1-fluid.res_gas)/(1-fluid.res_water)-...
            (sG)*fluid.dis_max*(fluid.res_water)/(1-fluid.res_water) )/fluid.dis_max;
        assert(all(rsH>-1e-4))
        rsH=max(rsH,0)
        %}
     else        
        smax=state.smax(:,2);
        rsH=state.s(:,1)*0;
     end
    diff=max(zt)-min(zt); 
    %subplot(2,2,1),cla
    %title('pressure')
    %plotCellData(G,state.pressure/barsa);colorbar('horiz'), caxis([100 600])
    
    %subplot(2,2,2),cla
    %title('saturation')
    %plotCellData(G,state.s(:,2));colorbar('horiz'), caxis([0 1])
    
    
    %%
    %zt=zt*0;
    figure(nnn)
     clf
     
    %
    
    if(true)
        p=state.pressure;
        sG=state.s(:,2);
        sGmax=state.sGmax;
        drho=fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p);
        h=(fluid.pcWG(sG,p,'sGmax',sGmax))./(drho*norm(gravity()));
        h_max=(fluid.pcWG(sGmax,p,'sGmax',sGmax))./(drho*norm(gravity()));
        %h_res1=(Gt.cells.H-(1-state.rs/fluid.dis_max).*Gt.cells.H);
        h_res=h_max+((state.rs.*state.s(:,1)-min_rs)./fluid.dis_max).*Gt.cells.H;
        %assert(all(abs((h_res1-h_res))<1e-6))
    else
        sG_free=free_sg(state.s(:,2),smax,opt);
        h=(sG_free.*Gt.cells.H)./(1-fluid.res_water);
        h_max=(sGmax.*Gt.cells.H)./(1-fluid.res_water);
        h_res=h_max+(-((sG.*fluid.res_water/(1-fluid.res_water)+(sGmax-sG)*(1-fluid.res_gas)/(1-fluid.res_water))*fluid.dis_max-state.s(:,1).*state.rs)).*Gt.cells.H./fluid.dis_max;
    end
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)], myCOColor(2))
  %
  %%{
    %sGmax=state.sGmax;
   
    
   
    
    %assert(all(h_res>h_max))
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_max(end:-1:1)], myCOColor(3))
    patch(xc([1:end end:-1:1]), ...
      [zt + h_max; zt(end:-1:1)+h_res(end:-1:1)], myCOColor(4))
    patch(xc([1:end end:-1:1]), ...
      [zt + h_res; zb(end:-1:1)], myCOColor(5))
  %} 
  assert(all(h_max>=h));
    %line(xc([1:end]),zt + h_max,'LineWidth',2,'Color','k')
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20], 0.7*[1 1 1]);
    set(gca,'YDir','reverse'), axis tight
    %%
    %axis([0 30 0 50])
    %print('-depsc2',[fname,'_fig_',num2str(nn),'panel.eps'])
    drawnow;   
    %pause(0.1)
%end
hold on
%axes(h1)

%
%{
clf,hold on
cartDims=Gt.cartDims;
nx=cartDims(1);
ind=1:nx+1;
x=Gt.nodes.coords(ind,1)/1e3;
z=Gt.nodes.z(ind);
x=[x;x(end:-1:1)];
z=[z;z(end:-1:1)+50];
patch(x,z,'y');set(gca,'YDir','reverse')
%set(gca,'FontSize',19)
%set(gca,'Ydir','reverse')
%}
set(gca,'FontSize',19)
set(gca,'Ydir','reverse')
xlabel('X km','FontSize',16)
ylabel('Depth m','FontSize',16)
set(gca,'Color','none')
%set(gca,'FontSize',16)

 %%

%rhoG =@(p) proj*fluid.rhoGS.*fluid.rhoGS(p*porl);
%rhoW =@(p) proj*fluid.rhoWS.*fluid.rhoWS(p*porl);
%p=state.pressure;
%sG=state.s(:,2);
%sGmax=state.sGmax;
%drho=fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p);
%h=(fluid.pcWG(sG,p,'sGmax',sGmax))./(drho*norm(gravity()));
%h_max=(fluid.pcWG(sGmax,p,'sGmax',sGmax))./(drho*norm(gravity()));
%h_res_b1=state.s(:,1).*(1-state.rs./fluid.dis_max).*Gt.cells.H;
%h_res_b=((state.rs.*state.s(:,1)-min_rs)./fluid.dis_max).*Gt.cells.H;
%plassert(all(abs(h_max+h_res_b+h_res_b1-Gt.cells.H)<1e-2))
%return
%%
h1 = gca;
if(nn==50)
    cells=[89,135,200,270]
elseif(nn==80)   
   cells=[130,270,350,405];  
else
    error()
   %cells=[200,270,300,400]; 
end
ha={};
for i=1:numel(cells)
cell=cells(i);
    pos= get(h1,'Position')
    ha{i}=axes()
if(i==1)   
    set(ha{i} ,'Position',get(h1,'Position')*0.3+[0.15,0.635,0.0, 0.0]);    
elseif(i==2)   
    set(ha{i} ,'Position',get(h1,'Position')*0.3+[0.30,0.13,0.0, 0.0]);        
elseif(i==3)   
    set(ha{i} ,'Position',get(h1,'Position')*0.3+[0.60,0.13,0.0, 0.0]);
elseif(i==4)    
    set(ha{i},'Position',get(h1,'Position')*0.3+[0.60,0.44,0.00, 0])    
else
    error()
    
end

%cell=100;
proj=zeros(Gt.cells.num,1);
proj(cell)=1;
prol=ones(Gt.cells.num,1);
%
%sco_max=0.7;
%%
%figure()
hold on
kscale=sqrt(rock.poro(cell)./(rock.perm(cell)))*fluidADI.surface_tension;
plot_hyst=true;
if(plot_hyst)
    h_int=h(cell)
    h_int_max=h_max(cell);
    h_int_res_b=h_res(cell);%h_int_max+h_res_b;%  Gt.cells.H(cell)-h_res_b(cell);
 else
    h_int=sco_max*H;
    h_int_max=sco_max*H; 
 end
 leg={};
 ss=[];hh=[];krk=[];ss_max=[];hh_max=[];
 
    [s,pc,kr,SH,krH, s_max,fval]= veRelpermTesterCell(h_int, drho(cell), fluid, H, 'samples',100,'hs_max',h_int_max,'kscale',kscale)
    ss=[ss,1-fval.s_h];
    ss_max=[ss_max,1-fval.s_hmax];
    hh_max=[hh_max,h_int_max-fval.h_max];
    krk=[krk,fval.kr_h];
    hh=[hh,h_int-fval.h];
    %plot(1-fval.s_h,h_int-fval.h,'LineWidth',2)
    %set(gca,'YDir','reverse')
    %axis([0 1,0 50])
 
 
 %plot(ss,hh,'LineWidth',2)
 %if(plot_hyst)
 %   plot(ss_max,hh_max,'--','LineWidth',2)
 %end
 s_max=@(d) interp1(hh_max,ss_max,d,'linear',0)
 if(h_int>0)
    s=@(d) interp1(hh,ss,d,'linear',0)
 else
    s=@(d) d*0; 
 end
 res_s =@(d) double(d<(h_int_res_b))
 w_s =@(d) ones(size(d));
 d=linspace(0,50,1345)';
 %
 ptch =@(x,y) [x;y(end:-1:1);x(1)];
 %ptch_h =@(x) [x;x(end:-1:1);x(1)];
 axes(ha{i})
 patch(ptch(s(d),s(d)*0),ptch(d,d),myCOColor(2))
 patch(ptch(s_max(d),s(d)),ptch(d,d),myCOColor(3))
 fix=max(res_s(d),s_max(d));
 patch(ptch(fix,s_max(d)),ptch(d,d),myCOColor(4))
 patch(ptch(w_s(d),fix),ptch(d,d),myCOColor(5))
if(~strcmp(surf_topo,'smooth'))
    line([0.0 1.0],[h_trap(cell) h_trap(cell)],'Color',myCOColor(1),'LineWidth',1)
end
 %legend(leg{:},'Location','NorthWest')
 set(gca,'YDir','reverse')
 %axis([0 1.02,0 4])
 box on
 set(gca,'LineWidth',2)
 set(gca,'FontSize',16)
 axes(h1)
 line([xc(cell),xc(cell)],[zt(cell)-50,zb(cell)+50],'LineWidth',4,'Color','r')
 axes(ha{i})
 %axx={h3}%,h4,h5}
%axx={h5}
end
top =@(x) interp1(xc,zt,x);
for i=1:numel(ha)
    axs=ha{i};
    %ax=axis(axs);
    %bpos=[ax([2,4]),ax([2,4])-ax([1,3])]
    %xpos=mean(ax([1,2])-0.03);
    xpos= xc(cells(i));
    bpos=[xpos,top(xpos)+25, 1 1];
    %pos5= get(h5,'Position');
    pos5=dsxy2figxy_new(axs,[0 25 1 1]);
    pos1=dsxy2figxy_new(h1,bpos +[0 0 0 0]);
    annotation('arrow',[pos1(1) pos5(1)],[pos1(2) pos5(2)],'LineWidth',2)
    axes(ha{i})
    axis off;
end
set(gcf,'PaperPositionMode','auto')
print('-depsc2',['figs/',fname,'fluid_',num2str(mf),'_',num2str(nn),'_with_detail.eps']);
end
end
end
return
plot(xc,(state.rs.*state.s(:,1)-min_rs)./fluid.dis_max,xc,(state.rs.*state.s(:,1)-(state.s(:,1)-(1-h_max./Gt.cells.H)).*fluid.dis_max)./fluid.dis_max)
% plot(xc,h_max+Gt.cells.H.*(state.rs.*state.s(:,1)-min_rs)./fluid.dis_max,xc,h_max+Gt.cells.H.*(state.rs.*state.s(:,1)-(state.s(:,1)-(1-h_max./Gt.cells.H)).*fluid.dis_max
 plot(xc,h_max./Gt.cells.H+(state.rs-min_rs./state.s(:,1))./fluid.dis_max,xc,1-(1-state.rs/fluid.dis_max).*state.s(:,1))