do_print=true;
if(do_print)
 mkdir('figs')
end
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
n_c=2;% plot top surface case
dissolution_c=true;%plot disolution case

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
zt=aquifer.Gt.cells.z;
zb=aquifer.Gt.cells.z+aquifer.Gt.cells.H;
Gt=aquifer.Gt;
for nn=[50,80]
    figure(nn),clf,hold on
    fprintf('Year from start injection %i year\n', floor(sum(res.control.step.val(1:nn-1))/year))
    fprintf('Year after stop of injection %i year\n', floor(sum(res.control.step.val(1:nn-1))/year)-50)
    state = res.results{kk}.states{nn};
    fluid = makeFluidModel(aquifer, 'residual', residual, ...
        'dissolution', dissolution, 'fluidType', 'sharp interface','only_pvt',false);%,'co2_type','coolprops');
    if(~dissolution)
        sG_free = free_sg(state.s(:,2),state.smax(:,2), ...
            struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
        sGmax=state.smax(:,2);
    else
        sG_free = free_sg(state.s(:,2),state.sGmax, ...
            struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
        sGmax=state.sGmax;
    end
    
    %%
    p=state.pressure;sG=state.s(:,2);
    rs=state.rs;   
    diff=max(zt)-min(zt); 
    
    
    h=(sG_free.*Gt.cells.H)./(1-fluid.res_water);
    
    patch(xc([1 1:end end]), [zt(end)-20; zt; zt(end)-20],.7*[1 1 1]);
    patch(xc([1 1:end end]), [zb(1)+20; zb; zb(1)+20],.7*[1 1 1]);
    patch(xc([1:end end:-1:1]), ...
      [zt + h_trap; zt(end:-1:1)], myCOColor(1))
    patch(xc([1:end end:-1:1]), ...
      [zt + h; zt(end:-1:1)+h_trap(end:-1:1)], myCOColor(2))
 
    h_max=(sGmax.*Gt.cells.H)./(1-fluid.res_water);
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
    
    % add innlet
    set(gca,'FontSize',19)
    set(gca,'Ydir','reverse')
    xlabel('X km','FontSize',16)
    ylabel('Depth m','FontSize',16)
    set(gca,'Color','none')
    h1 = gca;
    pos= get(h1,'Position');
    %delete(h3)
    h3=copyobj(h1,gcf);
    
    set(h3,'Position',get(h1,'Position')*0.3+[0.55,0.39,0.07, 0])
    axes(h3);%#ok
    set(gca,'FontSize',16)
    set(gca,'YDir','reverse');
    top =@(x) interp1(xc,zt,x);
    if(nn==50)
        xx=[10,10.5]+2.5;
    elseif(nn==80)
        xx=[21 23];
    else
        xx=[21 23];
        %error()
    end
    matop=max(top(xx));
    mitop=min(top(xx));
    if(nn==50)
        %matop=matop+50;
    end
    ax=[xx,mitop,matop];
    axis(ax);
    axes(h1);%#ok
    plot([ax(1),ax(2),ax(2),ax(1),ax(1)],[ax(3),ax(3),ax(4),ax(4),ax(3)],'r')
    axes(h3);%#ok
    axis off
    
    
    xlabel(''),ylabel('')
    h4=copyobj(h1,gcf);
    set(h4 ,'Position',get(h1,'Position')*0.3+[0.12,0.635,0.1, 0.0]);
    axes(h4)%#ok
    set(gca,'YDir','reverse');
    if(nn==50)
        xx=[5,5.5];
    elseif(nn==80)
        xx=[2.5 5];
    else
        xx=[2.5 5];
%        error()
    end
    
    matop=max(top(xx));
    if(nn>=80)
        matop=max(top(xx))+50;
    end
    mitop=min(top(xx));
    ax=[xx,mitop,matop];
    axis(ax);
    axes(h1);%#ok
    plot([ax(1),ax(2),ax(2),ax(1),ax(1)],[ax(3),ax(3),ax(4),ax(4),ax(3)],'r')
    axes(h3)%#ok
    axes(h4)%#ok
    axis off
    %%
    
    h5=copyobj(h1,gcf);
    %set(h4 ,'Position',get(h1,'Position')*0.35+[0.13,0.56,0.0, 0.0]);
    set(h5 ,'Position',get(h1,'Position')*0.3+[0.30,0.11,0.3, 0.0]);
    axes(h5)%#ok
    set(gca,'YDir','reverse');
    if(nn==50)
        xx=[7,8]+2;
    elseif(nn==80)       
        xx=[7.03,7.33]+9.3;        
    else
        xx=[7.03,7.33]+9.3;        
        %error()
    end
    
    matop=max(top(xx));
    if(nn==80)
        matop=max(top(xx));%+50;
    end
    mitop=min(top(xx));
    ax=[xx,mitop,matop];
    axis(ax);
    axes(h1);%#ok
    plot([ax(1),ax(2),ax(2),ax(1),ax(1)],[ax(3),ax(3),ax(4),ax(4),ax(3)],'r')
    axes(h3);%#ok
    axes(h4);%#ok
    axes(h5);%#ok
    axis off
    %%
    axx={h3,h4,h5};
    for i=1:numel(axx)
        axs=axx{i};
        ax=axis(axs);
        bpos=[ax([2,4]),ax([2,4])-ax([1,3])];
        xpos=mean(ax([1,2])-0.03);
        bpos([1,2])=[xpos,top(xpos)];
        %pos5= get(h5,'Position');
        pos5=dsxy2figxy_new(axs,bpos -[0 0 0 0]);
        pos1=dsxy2figxy_new(h1,bpos +[0 0 0 0]);
        annotation('arrow',[pos1(1) pos5(1)],[pos1(2) pos5(2)],'LineWidth',2)
    end   
    %%
    set(gcf,'PaperPositionMode','auto')   
    %axis([0 25 0 50])
    if(do_print)        
        print('-depsc2',['figs/fig14_',num2str(nn),'_detail.eps'])
    end
    drawnow;   
    pause(0.1)
end

    
    
    

