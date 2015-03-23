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
kk=1;

 for n=1:2, % flat or non flat topsurface
   %for residual= [false,true] %residual saturation or not
    residual=true;
    %
    f=figure(),clf,hold on;
    set(f,'Position',[0,600,700,500]);
    k=1; 
    %ivec=[1,2,3,4];
     for i= 1:4 %residual saturation or not   
        %% Make fluid model
        for dissolution=[false,true]
            state = res.results{kk}.states{end-70};
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', fluid_types{i},'only_pvt',true);%,'co2_type','coolprops');
            if(~dissolution)
                sG = free_sg(state.s(:,2),state.smax(:,2), ...
                    struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
            else
                sG = free_sg(state.s(:,2),state.sGmax, ...
                    struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
            end            
            hold on
            if(dissolution)
               mline=[line_colors{i},'--'] ;
            else
              mline=[line_colors{i},'-'] ;
            end
            ff=res.results{9}.ff;            
            plot(res.xc, filter2(ff,sG.*Gt.cells.H), mline, 'LineWidth', 2);    
            %fff=1;plot(res.xc, filter2(fff,sG.*Gt.cells.H), mline, 'LineWidth', 1);    
            fprintf('Residual %f %f\t',fluid.res_gas, fluid.res_water)
            fprintf(' flat %i \t', n)
            if(dissolution)
                fprintf('disolution %f %f\n',fluid.dis_max, fluid.dis_rate)
            else
                fprintf('\n');
            end
            %hold off
            drawnow;                        
            k = k+1;
            kk=kk+1;
        end
      
     end
    if(n==2)
        plot(res.xc(50:end-50),hts(50:end-50),'r-','LineWidth',2)
    end
    
    axis tight
    ax=axis();
    %axis([ 0 25 ax(3) ax(4)])
    set(gca,'YDir','reverse','FontSize',16);
    % Add to legends one for color and one for linestyle
    %axis([0 5 0 1])
    ivec=[1,2,4,3];
    if(n==1)
        addLegends(gca, {fluid_types{ivec}}, {line_colors{ivec}},{legendtext{1:2}},{'-','--'});
    else
        addLegends(gca, {fluid_types{ivec}}, {line_colors{ivec}},{legendtext{3:4}},{'-','--'});
    end
    %legend(legendtext{:}, 4);
    if(do_print)
        if(n==1)
            print -depsc2 figs/ex1-fig21a.eps;
        else
            print -depsc2 figs/ex1-fig21b.eps;
        end
    end
end
