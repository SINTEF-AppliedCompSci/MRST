%%
smooth=false;
depth=1300;
res_fluid=true;
mydir='test_data_dis_all/';fluid.dis_max=0.03;
mydir='data_all_dis/';
mydir='data/';opt.dis_max=0.03;
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
surf_topo='inf rough';example1D_makeparam; % only to get htrap
Gt = res.opt.Gt;


%%
%% Plot results
%figure
states=res.results{1}.states;
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
%fluid.res_gas=sr;
%fluid.res_oil=sw;
opt.res_gas=sr;
opt.res_oil=sw;

%if(use_dis)
%   fluid.dis_max=0.03;
%end
%fluid=fluid{1};

%for nn=1:numel(states)
for nn=[50,80]
    nnn=nn
    state=states{nn};
    p=state.pressure;sG=state.s(:,2);sGmax=state.sGmax;
    rs=state.rs;
    if(isfield(fluid,'dis_max'))
        smax=state.sGmax;
        sGmax=state.sGmax;
        sG=state.s(:,2);
        %rsH=Gt.cells.H.*state.s(:,1).*state.rs/fluid.dis_max;
        rsH=Gt.cells.H.*(state.s(:,1).*state.rs-...
            (sGmax-sG)*fluid.dis_max*(1-fluid.res_gas)/(1-fluid.res_oil)-...
            (sG)*fluid.dis_max*(fluid.res_oil)/(1-fluid.res_oil) )/fluid.dis_max;
        assert(all(rsH>-1e-4))
        rsH=max(rsH,0)
    else
        smax=state.smax(:,2);
        rsH=state.s(:,1)*0;
    end
    diff=max(zt)-min(zt);
    
    
    %filter2(ff./sum(ff),hh))
    %%
    %zt=zt*0;
    figure(nnn)
    clf
    sG=free_sg(state.s(:,2),smax,opt);
    %
    h=(sG.*Gt.cells.H)./(1-fluid.res_oil);
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
        [zt + h_trap; zt(end:-1:1)], myCOColor(1))
    patch(xc([1:end end:-1:1]), ...
        [zt + h; zt(end:-1:1)+h_trap(end:-1:1)], myCOColor(2))
    %
    %%{
    %sGmax=state.sGmax;
    
    h_max=(sGmax.*Gt.cells.H)./(1-fluid.res_oil);
    h_res=(-((sG.*opt.res_oil/(1-fluid.res_oil)+(sGmax-sG)*(1-opt.res_gas)/(1-opt.res_oil))*opt.dis_max-state.s(:,1).*state.rs)).*Gt.cells.H./opt.dis_max;
    assert(all(h_res>-1e-4))
    patch(xc([1:end end:-1:1]), ...
        [zt + h; zt(end:-1:1)+h_max(end:-1:1)],myCOColor(3))
    patch(xc([1:end end:-1:1]), ...
        [zt + h_max; zt(end:-1:1)+h_max(end:-1:1)+h_res(end:-1:1)], myCOColor(4))
    patch(xc([1:end end:-1:1]), ...
        [zt + h_max+h_res; zb(end:-1:1)], myCOColor(5))
    %}
    line(xc([1:end]),zt + h_max,'LineWidth',2,'Color','k')
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    set(gca,'YDir','reverse'), axis tight
    %%
    %axis([0 30 0 50])
    drawnow;
    pause(0.1)
    %end
    hold on
    
    set(gca,'FontSize',19)
    set(gca,'Ydir','reverse')
    xlabel('X km','FontSize',16)
    ylabel('Depth m','FontSize',16)
    set(gca,'Color','none')
    h1 = gca;
    pos= get(h1,'Position')
    %delete(h3)
    h3=copyobj(h1,gcf);
    
    set(h3,'Position',get(h1,'Position')*0.3+[0.55,0.39,0.07, 0])
    axes(h3)
    set(gca,'FontSize',16)
    set(gca,'YDir','reverse');
    top =@(x) interp1(xc,zt,x);
    if(nn==50)
        xx=[5,6];
        xx=[10,10.5]+2.5;
    elseif(nn==80)
        xx=[21 23];
    else
        error()
    end
    matop=max(top(xx));
    mitop=min(top(xx));
    if(nn==50)
        %matop=matop+50;
    end
    ax=[xx,mitop,matop];
    axis(ax);
    axes(h1)
    plot([ax(1),ax(2),ax(2),ax(1),ax(1)],[ax(3),ax(3),ax(4),ax(4),ax(3)],'r')
    axes(h3)
    axis off
    
    
    xlabel(''),ylabel('')
    h4=copyobj(h1,gcf);
    set(h4 ,'Position',get(h1,'Position')*0.3+[0.12,0.635,0.1, 0.0]);
    axes(h4)
    set(gca,'YDir','reverse');
    if(nn==50)
        xx=[5,5.5];
    elseif(nn==80)
        xx=[2.5 5];
    else
        error()
    end
    
    matop=max(top(xx));
    if(nn==80)
        matop=max(top(xx))+50;
    end
    mitop=min(top(xx));
    ax=[xx,mitop,matop];
    axis(ax);
    axes(h1)
    plot([ax(1),ax(2),ax(2),ax(1),ax(1)],[ax(3),ax(3),ax(4),ax(4),ax(3)],'r')
    axes(h3)
    axes(h4)
    axis off
    %%
    
    h5=copyobj(h1,gcf);
    %set(h4 ,'Position',get(h1,'Position')*0.35+[0.13,0.56,0.0, 0.0]);
    set(h5 ,'Position',get(h1,'Position')*0.3+[0.30,0.11,0.3, 0.0]);
    axes(h5)
    set(gca,'YDir','reverse');
    if(nn==50)
        xx=[7,8]+2;
    elseif(nn==80)
        xx=[7,7.2]+4;% non flowing
        xx=[7.03,7.33]+9.3;
        
    else
        error()
    end
    
    matop=max(top(xx));
    if(nn==80)
        matop=max(top(xx));%+50;
    end
    mitop=min(top(xx));
    ax=[xx,mitop,matop];
    axis(ax);
    axes(h1)
    plot([ax(1),ax(2),ax(2),ax(1),ax(1)],[ax(3),ax(3),ax(4),ax(4),ax(3)],'r')
    axes(h3)
    axes(h4)
    axes(h5)
    axis off
    %%
    axx={h3,h4,h5}
    %axx={h5}
    for i=1:numel(axx)
        axs=axx{i};
        ax=axis(axs);
        bpos=[ax([2,4]),ax([2,4])-ax([1,3])]
        xpos=mean(ax([1,2])-0.03);
        bpos([1,2])=[xpos,top(xpos)]
        %pos5= get(h5,'Position');
        pos5=dsxy2figxy_new(axs,bpos -[0 0 0 0]);
        pos1=dsxy2figxy_new(h1,bpos +[0 0 0 0]);
        annotation('arrow',[pos1(1) pos5(1)],[pos1(2) pos5(2)],'LineWidth',2)
    end
    %axes(h3)
    %axes(h4)
    %axes(h5)
    %%
    set(gcf,'PaperPositionMode','auto')
    print('-depsc2',['figs/',fname,'_',num2str(nn),'with_detail.eps']);
end
