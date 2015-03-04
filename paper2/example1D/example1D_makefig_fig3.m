%%
ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
ff_names={ff_names{2:end}}
ff_names={'fine scale','square','inf rough','sinus'}
results={};
res={};
fnames={};cases={};
mydir='data_smotted_all_2/';% include only runes with linear relperm
mydir='data/'
%mydir='data_example_1D/';
%{
for smooth=[true,false]
for res_fluid=[false,true]
for depth=[2300,1300]
%}
%%{% the linear with res fluid seems to be a bit unstable
%for smooth=[true]
surf_topos={'smooth','square','inf_rough','sinus'}
for k=[1,2,3,numel(surf_topos)]
    surf_topo=surf_topos{k};   
for res_fluid=[true]
for depth=[2300]
     if(res_fluid)
            fname=['topo_model_',surf_topo,'_res_fluid','_','_depth_',num2str(depth)];
        else
            fname=['topo_model_',surf_topo,'_nores_fluid','_','_depth_',num2str(depth)];
     end
   res{end+1}=load([mydir,fname,'.mat']);
   sims={res{end}.results{:}};
   for dd=1:numel(sims)
    cases{end+1}=struct('fname',fname,'depth',depth,'res_fluid',res_fluid,'surf_topo',surf_topo);
   end
   %results={results{:},res{end}.results{:}};
   results={results{:},sims{:}};
end
end
end

%%

xc=res{end}.xc
xx=xc(150:650);ff=exp(-((xx-xc(400))/(0.3)).^2);plot(xx,ff)
ff=ff/sum(ff);
ff=1;
s=results{1}.states{end}.s(:,2)
clf,plot(xc,s,xc,filter2(ff./sum(ff),s))
%%
figure(1),clf,hold on
%
hold on;
marker={'s','*','d','o','<'}
marker={'ks','b*','md','ro','<'}
markers={'k','g','m','r','b'}
marker={};
for i=1:numel (results)
    marker={marker{:},markers{:}}
end
use_filter=true;
zt = res{end}.z;
zb=zt+50;H=50;
for k=1:numel(results)
    if(strcmp(cases{k}.surf_topo,'smooth'))
        xx=xc(150:650);ff=exp(-((xx-xc(400))/(0.3)).^2);%plot(xx,ff)
        ff=ff/sum(ff);
        %ff=1;
    else
        ff=1;
    end
    
    
    states=results{k}.states;
    %p0=x0.pressure;
    %xc = Gt.cells.centroids(:,1)/1e3;
    
    %zb = zt + Gt.cells.H;
    ns=numel(states);
    %for nn=[numel(dt),ns-10]%ns];%[numel(dt_inj),ns]
    sss=70;
    for nn=[ns-sss:ns-sss]%ns];%[numel(dt_inj),ns]    
    %for nn=[ns:ns]%ns];%[numel(dt_inj),ns]        
    %for nn=[ns-100:20:ns-10]%ns];%[numel(dt_inj),ns]        
        %for nn=ns:ns
        
        state=states{nn};
        
        if(false )%isfield(fluid,'dis_max'))
            smax=state.sGmax;
            rsH=state.s(:,1).*state.rs/fluid.dis_max;
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
        
        % plot as VE
        %{
    subplot(3,1,2),hold on
    plot(xc,state.pressure/barsa); set(gca,'YLim',[100 300]);
    subplot(3,1,3),hold on
    plot(xc,(state.pressure-p0)/barsa); %set(gca,'YLim',[0 200]);
    subplot(3,1,1),hold on
        %}
        if(cases{k}.res_fluid)
            %sr= 0.21;sw= 0.11
           opt.res_gas=0.21;
           opt.res_oil=0.11; 
        else
           opt.res_gas=0;
           opt.res_oil=0;
        end
        sG=free_sg(state.s(:,2),smax,opt);
        if(use_filter==false)
            if(isfield(fluid,'dis_max'))
                plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H,xc,smax.*Gt.cells.H+rsH); %set(gca,'YLim',[0 20]);
            else
                %plot(xc,state.s(:,2).*Gt.cells.H,['-',marker{k}],xc,smax.*Gt.cells.H,['-',marker{k}]);
                plot(xc,sG.*H,['-',marker{k}],'LineWidth',2)
                plot(xc,state.s(:,2).*H,['.',marker{k}],'LineWidth',2)
                %plot(xc,smax.*Gt.cells.H,['-',marker{k}],'LineWidth',2) 
            end
        else
           plot(xc,filter2(ff,sG.*H),['-',marker{k}],'LineWidth',2)
           %plot(xc,filter2(ff,smax.*Gt.cells.H),['.',marker{k}],'LineWidth',2) 
           %plot(xc,(sG.*Gt.cells.H))%,['-',marker{k}],'LineWidth',2)
           %plot(xc,(smax.*Gt.cells.H))%,['.',marker{k}],'LineWidth',2) 
        end
        %%
    end
    %%
end
set(gca,'FontSize',19)
set(gca,'Ydir','reverse')
h1 = gca;
h3 = axes('Position',get(h1,'Position'));
x=1:2;
leg_names=ff_names;
tmp=nan(numel(x),numel(leg_names));
%tmp=repmat(x,numel(leg2),1);
f2=plot(x,tmp)
axis off;
set(gca,'XtickLabel',[])
for f=1:numel(f2)
    set(f2(f),'Marker','none');
    set(f2(f),'MarkerFaceColor','none');
    set(f2(f),'LineStyle','-');
    %set(f2(f),'Color',get(a(f),'Color'));
    %set(f2(f),'Color',mycol(f,:));
    set(f2(f),'Color',marker{f});
end
set(h3,'FontSize',16);
bb=legend(leg_names,'Location','SouthEast','FontSize',16);
set(gca,'FontSize',19)

print('-depsc2','figs/example1D_fig3.eps')
return
