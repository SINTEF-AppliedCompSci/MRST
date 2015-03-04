%%
close all
%ff_names={'sharp interface','linear cap.','P-scaled table','P-K-scaled table','S table'};
ff_names={'no res','res'};
%ff_names={ff_names{2:end}}
results={};
res={};
fnames={};cases={};
mydir='data_example_1D_linear/';% include only runes with linear relperm
mydir='data/'
for smooth=[true]%,false]
for res_fluid=[false,true]
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
   sims={res{end}.results{1}};
   for dd=1:numel(sims)
    cases{end+1}=struct('fname',fname,'depth',depth,'smooth',smooth,'res_fluid',res_fluid);
   end
   results={results{:},sims{:}};
end
end
end

%%

%
hold on;
markers={'k','r'}

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
        %ff=1;
    else
        ff=1;
    end
    
    
    states=results{k}.states; 
    ns=numel(states);   
    sss=70
    for nn=[ns-sss,ns-sss]
        state=states{nn};
        if(false)
            smax=state.sGmax;
            rsH=state.s(:,1).*state.rs/fluid.dis_max;
        else
            smax=state.smax(:,2);
            rsH=state.s(:,1)*0;
        end
        diff=max(zt)-min(zt);
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
                plot(xc,sG.*Gt.cells.H,['-',marker{k}],'LineWidth',2)
                %plot(xc,state.s(:,2).*Gt.cells.H,['.',marker{k}],'LineWidth',2)
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
print('-depsc2','figs/example1D_fig1')

return

