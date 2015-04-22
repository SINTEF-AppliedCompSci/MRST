%%
close all
ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
ff_names={ff_names{1:end}}
ff_names={'no res. flat','res. flat','no res. oss.','res. oss.'}
results={};
res={};
fnames={};cases={};
mydir='data_example_1D_linear/';% include only runes with linear relperm
mydir='data/'
for smooth=[false]
for res_fluid=[false]%,true]
for depth=[2300]
%}
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
   res{end+1}=load([mydir,fname,'.mat']);
   %sims={res{end}.results{2:end}};
   sims={res{end}.results{1}};
   for dd=1:numel(sims)
    cases{end+1}=struct('fname',fname,'depth',depth,'smooth',smooth,'res_fluid',res_fluid);
   end
   %results={results{:},res{end}.results{:}};
   results={results{:},sims{:}};
end
end
end

%%


figure()
hold on
markers={'--b','b','--r','r'}
markers={'k','k','b','r'}
marker={};
for i=1:numel (results)
    marker={marker{:},markers{:}}
end
use_filter=true;
zt = res{end}.z;
zb=zt+50;H=50;
for k=1:numel(results)
    if(~cases{k}.smooth)
        xx=xc(150:650);ff=exp(-((xx-xc(400))/(0.3)).^2);%plot(xx,ff)
        ff=ff/sum(ff);
        ff=1;
    else
        ff=1;
    end
    
    
    states=results{k}.states;
    %p0=x0.pressure;
    %xc = Gt.cells.centroids(:,1)/1e3;
    %zt = Gt.cells.z;
    %zb = zt + Gt.cells.H;
    ns=numel(states);
    %for nn=[numel(dt),ns-10]%ns];%[numel(dt_inj),ns]
    sss=70
    for nn=[ns-sss,ns-sss]%ns];%[numel(dt_inj),ns]    
    %for nn=[ns:ns]%ns];%[numel(dt_inj),ns]        
    %for nn=[ns-100:20:ns-10]%ns];%[numel(dt_inj),ns]        
        %for nn=ns:ns
        
        state=states{nn};
        
        if(false)%isfield(fluid,'dis_max'))
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
           opt.res_water=0.11; 
        else
           opt.res_gas=0;
           opt.res_water=0;
        end
        sG=free_sg(state.s(:,2),smax,opt);
        if(use_filter==false)
            if(isfield(fluid,'dis_max'))
                plot(xc,state.s(:,2).*Gt.cells.H,xc,smax.*Gt.cells.H,xc,smax.*Gt.cells.H+rsH); %set(gca,'YLim',[0 20]);
            else
                %plot(xc,state.s(:,2).*Gt.cells.H,['-',marker{k}],xc,smax.*Gt.cells.H,['-',marker{k}]);
                plot(xc,sG.*Gt.cells.H,['-',marker{k}],'LineWidth',2)
                %plot(xc,state.s(:,2).*Gt.cells.H,['.',marker{k}],'LineWidth',2)
                %plot(xc,smax.*Gt.cells.H,['-',marker{k}],'LineWidth',2) 
            end
        else
           %aa(k)=
            plot(xc,filter2(ff,sG.*H),[marker{k}],'LineWidth',1)
           %plot(xc,filter2(ff,smax.*Gt.cells.H),['.',marker{k}],'LineWidth',2) 
           %plot(xc,(sG.*Gt.cells.H))%,['-',marker{k}],'LineWidth',2)
           %plot(xc,(smax.*Gt.cells.H))%,['.',marker{k}],'LineWidth',2) 
        end
        %%
    end
    %%
end
%%
clf
x=xc(150:650);ff=exp(-((xx-xc(400))/(0.3)).^2);%plot(xx,ff)
ff=ff/sum(ff);
plot(xc,filter2(ff,sG.*H),'k','LineWidth',2),hold on
plot(xc,sG.*H,'k','LineWidth',1)
set(gca,'YDir','reverse')
%%
h1 = gca;
pos= get(h1,'Position')
h3 = axes('Position',get(h1,'Position')*0.3+[0.6,0.6,0, 0]);
plot(xc,filter2(ff,sG.*H),'k','LineWidth',2),hold on
plot(xc,sG.*H,'k','LineWidth',1)
set(gca,'YDir','reverse');axis([15 16 0 6])
set(gca,'XTick',[15 15.5 16])
axes(h1)
plot([15 16 16 15 15]',[0.5 0.5 5 5 0.5],'r','LineWidth',2)
axes(h3)
%%
x=res{end}.xc;
z=res{end}.z;
x=[x;x(end:-1:1)];
z=[z;z(end:-1:1)+50];

h4 = axes('Position',get(h1,'Position')*0.3+[0.6,0.2,0, 0]);
patch(x,z,'b');set(gca,'YDir','reverse');
x=res{end}.xc;
z=res{end}.z;
x=[x;x(end:-1:1)];
z=[z;z(end:-1:1)+sG(end:-1:1).*H];
patch(x,z,'r')
axis([15 16 1895 1925])
%axis([15 16 1880 1970])
%set(h4,'YAxisLocation','right') 
set(gca,'XTick',[15 15.5 16])
%set(h4,'XTick',[25 26])
set(gcf,'PaperPositionMode','auto')
%set(h3,'XTick',[5 6]);set(h3,'YTick',[1150 1175 1200])






%%



print('-depsc2','figs/example1D_fig2_b')
return
