do_print=true;
if(do_print)
 mkdir('figs')
end
gravity on;


depth=1300;
%res=load('data/disolutionExample1Data');
res=load(['data/disolutionExample1Data_',num2str(depth),'.mat']);
legendtext = {'No dissolution (A=0)', 'Dissolution (A=0)', ...
    'No dissolution (A=2)', 'Dissolution (A=2)'};
linetype = {'b-', 'b--', 'r-', 'r--'};
aquifer = makeAquiferModel('A',0,'D',depth);

G  = aquifer.G;
Gt = aquifer.Gt;
z  = G.cells.centroids(:,3);
% left and right panel of figure 15
%%
kk=1;
for residual= [false,true] %residual saturation or not   
    f=figure(),clf,
    set(f,'Position',[0,600,700,500]);
    k=1;
    for n=1:2, % flat or non flat topsurface
        %% Make fluid model
        for dissolution=[false,true]
            state = res.results{kk}.states{70};
            %state = res.results{kk}.states{end-80};
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', 'sharp interface','only_pvt',true);%,'co2_type','coolprops');  
            if(~dissolution)
                sG = free_sg(state.s(:,2),state.smax(:,2), ...
                    struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
            else
                sG = free_sg(state.s(:,2),state.sGmax, ...
                    struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
            end                                   
            hold on
            plot(res.xc, filter2(res.results{kk}.ff,sG.*Gt.columns.dz), linetype{k}, 'LineWidth', 2);
            fprintf('Residual %f %f\t',fluid.res_gas, fluid.res_water)
            fprintf(' flat %i \t', n)
            if(dissolution)
                fprintf('disolution %f %f\n',fluid.dis_max, fluid.dis_rate)
            else
                fprintf('\n');
            end
            hold off
            drawnow;                        
            k = k+1;
            kk=kk+1;
        end
    end
    axis tight
    ax=axis();
    axis([ 0 28 ax(3) ax(4)])
    set(gca,'YDir','reverse','FontSize',16);
   legend(legendtext{:}, 4,'Location','SouthWest');
    if(do_print)
        if(~residual)
            print -depsc2 figs/ex1-fig15a.eps;
        else
            print -depsc2 figs/ex1-fig15b.eps;
        end
    end
end
