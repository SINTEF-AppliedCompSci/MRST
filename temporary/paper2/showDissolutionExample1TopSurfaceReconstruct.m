do_print=true;
if(do_print)
    mkdir('figs')
end
gravity on;
res=load('data/disolutionTopSurfaceExample1Data.mat');
%res=b;
%res=a;
legendtext = {'No dissolution (A=0)', 'Dissolution (A=0)', ...
    'No dissolution (A=2)', 'Dissolution (A=2)'};
%linetype = {'b-', 'b--', 'r-', 'r--'};
line_colors={'k','b','m','g'};
upaquifer = makeAquiferModel('A',0);
faquifer = makeAquiferModel('A',0);

upAquifer = makeAquiferModel('A',0);
fAquifer  = makeAquiferModel('A',2);
aquifer=upAquifer;

z  = fAquifer.Gt.cells.z;
zt = max(z)*ones(size(z));
for i=2:numel(zt)-1
    zt(i)=max(z(i:end));
end
zt(end) = max(zt(end-1),z(end));
ht  = zt - z;
ff = exp(-linspace(-25,25,501).^2); ff=ff'/sum(ff);
hts = filter2(ff, ht);
%G  = aquifer.G;
Gt = aquifer.Gt;
z  = Gt.cells.z;
%z  = G.cells.centroids(:,3);
% left and right panel of figure 15
%%
fluid_types={'sharp interface','linear cap','P-scaled table','P-K-scaled table'};

dissolution_c=true;n_c=2;
for i_c=1:3;






residual_c=true;
k=1;
for n=1:2, % flat or non flat topsurface
    %for residual= [false,true] %residual saturation or not
    residual=residual_c;
        %ivec=[1,2,3,4];
        for i= 1:4 %residual saturation or not
            for dissolution=[false,true]
                if(n_c==n && i_c==i && residual_c==residual && dissolution_c==dissolution)
                    kk=k;
                end
                k=k+1;
            end
        end
    %end
end

if(n_c==2)
    aquifer=fAquifer;
    h_trap=ht;
else
    aquifer=upAquifer;
    h_trap=0*aquifer.Gt.cells.H;
end
Gt=aquifer.Gt;

for nn=[50,80]
state = res.results{kk}.states{nn};
fluid = makeFluidModel(aquifer, 'residual', residual_c, ...
    'dissolution', dissolution_c, 'fluidType', fluid_types{i_c},'only_pvt',false);%,'co2_type','coolprops');

%%{
if(~dissolution_c)
    sG_free = free_sg(state.s(:,2),state.smax(:,2), ...
        struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
else
    sG_free = free_sg(state.s(:,2),state.sGmax, ...
        struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
end
%}
% plot
f=figure(),clf,hold on;
set(f,'Position',[0,600,700,500]);
p=state.pressure;
sG=state.s(:,2);
if(dissolution_c)
    sGmax=state.sGmax;
else
    sGmax=state.smax(:,2);
end
drho=fluid.rhoWS.*fluid.bW(p)-fluid.rhoGS.*fluid.bG(p);
h=(fluid.pcWG(sG,p,'sGmax',sGmax))./(drho*norm(gravity()));
h_max=(fluid.pcWG(sGmax,p,'sGmax',sGmax))./(drho*norm(gravity()));
if(dissolution_c)
        mm=minRs(p,sG,sGmax,fluid,Gt);%.*Gt.cells.H;
        h_res_diff=Gt.cells.H.*((1-sG).*state.rs-mm)/fluid.dis_max;
    else
        h_res_diff=0*Gt.cells.H;
end
%h_res_diff=((state.rs.*state.s(:,1)-min_rs)./fluid.dis_max).*Gt.cells.H;
assert(all(h_res_diff>=-1e-5));
h_res=h_max+h_res_diff;

%
zt=aquifer.Gt.cells.z;
zb=zt+aquifer.Gt.cells.H;
xc=aquifer.Gt.cells.centroids(:,1)/1e3;
    
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    
    patch(xc([1:end end:-1:1]), ...
      [zt + h_trap; zt(end:-1:1)], myCOColor(1))
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_trap(end:-1:1)], myCOColor(2))
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_max(end:-1:1)], myCOColor(3))
    patch(xc([1:end end:-1:1]), ...
      [zt + h_max; zt(end:-1:1)+h_res(end:-1:1)], myCOColor(4))
    patch(xc([1:end end:-1:1]), ...
      [zt + h_res; zb(end:-1:1)], myCOColor(5))
      assert(all(h_max>=h));
    %line(xc([1:end]),zt + h_max,'LineWidth',2,'Color','k')
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20], 0.7*[1 1 1]);
    
    set(gca,'YDir','reverse'), axis tight

    drawnow;   
hold on



set(gca,'FontSize',19)
set(gca,'Ydir','reverse')
xlabel('X km','FontSize',16)
ylabel('Depth m','FontSize',16)
set(gca,'Color','none')
%set(gca,'FontSize',16)


h1 = gca;
if(nn==50)
    cells=[89,135,200,270,89]
elseif(nn==80)
    %cells=[129,270,350,405,129];
    cells=[89,270,350,405,89];
else
    cells=[89,135,200,270,89];%[130,270,350,405];
    %error()
    %cells=[200,270,300,400];
end
if(n_c==1)
    cells=cells(1:4);
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
   elseif(i==5)
        set(ha{i},'Position',get(h1,'Position')*0.3+[0.4,0.78,-0.05, -0.15])   
    else
        error()
    end
    
    proj=zeros(Gt.cells.num,1);
    proj(cell)=1;
    prol=ones(Gt.cells.num,1);
    %
    %sco_max=0.7;
    %%
    %figure()
    hold on
    rock=aquifer.rock2D;
    kscale=sqrt(rock.poro(cell)./(rock.perm(cell)))*fluid.surface_tension;
    plot_hyst=true;
    H=aquifer.Gt.cells.H(cell);
    if(residual_c)
        h_int=h(cell)
        h_int_max=h_max(cell);
        h_int_res_b=h_res(cell);%h_int_max+h_res_b;%  Gt.cells.H(cell)-h_res_b(cell);
    else
        h_int=h(cell);
        h_int_max=h(cell);
    end
    leg={};
    ss=[];hh=[];krk=[];ss_max=[];hh_max=[];
    
    [s,pc,kr, s_max, s_free, fval]= veRelpermTesterCell(h_int, drho(cell), fluid, H, 'samples',1000,'hs_max',h_int_max,'kscale',kscale);
    mtol=1e-2;
    assert(abs(s_free-sG_free(cell))<mtol);
    assert(abs(s-sG(cell))<mtol);
    assert(abs(s_max-sGmax(cell))<mtol);
    pc_cells=fluid.pcWG(sG,p,'sGmax',sGmax);
    assert(abs(pc-pc_cells(cell))<mtol);
    krG=fluid.krG(sG,p,'sGmax',sGmax);
    assert(abs(kr-krG(cell))<mtol);
    
    ss=[ss,1-fval.s_h];
    ss_max=[ss_max,1-fval.s_hmax];
    hh_max=[hh_max,h_int_max-fval.h_max];
    krk=[krk,fval.kr_h];
    hh=[hh,h_int-fval.h];
    
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
    if(false)
    patch(ptch(s(d),s(d)*0),ptch(d,d),myCOColor(2))
    patch(ptch(s_max(d),s(d)),ptch(d,d),myCOColor(3))
    fix=max(res_s(d),s_max(d));
    patch(ptch(fix,s_max(d)),ptch(d,d),myCOColor(4))
    patch(ptch(w_s(d),fix),ptch(d,d),myCOColor(5))
    if(n_c==2)
        line([0.0 1.0],[h_trap(cell) h_trap(cell)],'Color',myCOColor(1),'LineWidth',1)
    end
    
    else
        assert(all(s(d)>=0))
      s_res=(s_max(d)-s(d)).*fluid.res_gas;%'non flowing co2'  
      patch(ptch(s_res,s(d)*0),ptch(d,d),myCOColor(3))  
      s_real=s(d)+s_res;%total saturation
      ind=find(d>h_trap(cell),1,'first');
      patch(ptch(s_real(1:ind),s_res(1:ind)),ptch(d(1:ind),d(1:ind)),myCOColor(1))        
      patch(ptch(s_real(ind+1:end),s_res(ind+1:end)),ptch(d(ind+1:end),d(ind+1:end)),myCOColor(2))        
      
      %fix=max(res_s(d),s_max(d));
      fix=res_s(d);
      patch(ptch(fix,s_real),ptch(d,d),myCOColor(4))
      patch(ptch(w_s(d),fix),ptch(d,d),myCOColor(5))
      if(n_c==2)
            line([0.0 1.0],[h_trap(cell) h_trap(cell)],'Color',myCOColor(1),'LineWidth',1)
      end
      ss_max=s_max(d);
      ss_max(d>h_int_max)=nan;
      ind=find(d>h_int_max,1,'first')
      ss_max(ind)=0;
      plot(ss_max,d,'k');%,'LineWidth',2) 
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
    if(i<5)
        xpos= xc(cells(i));
        bpos=[xpos,top(xpos)+25, 1 1];
        %pos5= get(h5,'Position');
        pos5=dsxy2figxy_new(axs,[0 25 1 1]);
        pos1=dsxy2figxy_new(h1,bpos +[0 0 0 0]);
        annotation('arrow',[pos1(1) pos5(1)],[pos1(2) pos5(2)],'LineWidth',2)
    else
        if(n_c==2)
            axes(ha{i})
            ax=axis();
            axis([0 1,0,h_trap(cells(i))*5]);
        else
           %axes(ha{i})
           %ax=axis();
           %axis([0 1,0,h(cells(i))*5]); 
        end
    end
    axes(ha{i})
    axis off;
end
set(gcf,'PaperPositionMode','auto')
if(do_print)
    if(nn==50)
        print('-depsc2',['figs/ex1-fig22_row1_panel_',num2str(i_c),'.eps']) 
    else
      print('-depsc2',['figs/ex1-fig22_row2_panel_',num2str(i_c),'.eps'])   
    end
end
end
end
