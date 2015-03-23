%% Example x: Compare the effects small-scale undulations
% In this example we consider a 1D antiform aquifer with a caprock given by
% the following expression
%      z = D - L1 * sin(x/L1)*tan(phi) + A sin(2*pi*x/L2)
% Injection of CO2 is simulate using models with/without residual
% saturation. This is done for the case of a flat surface (A=0) and with
% the case with small-scale caprock undulations (A=2). For all this cases
% we simulate with and without disolution effects.  
% We show that the disolution of the residual saturation prevent the disolution 
% of the free CO2 left in the small traps in our model which first disolve the
% residual sone. Finaly we show the distribution of the CO2 in real space.
% The data from the simulations are also stored and can be inspected the
% figures can be generate by showDissolutionExample1.m
try
    require co2lab ad-fi
catch %#ok<CTCH>
    mrstModule add co2lab ad-fi
end
gravity reset on
do_print=true;mkdir('figs')
%% Time steps for injection and migration
Ti  =   50*year;
dTi =  2*year;
istep = linspace(0.1*year, dTi, 10)';
istep = [istep; ones(floor((Ti-sum(istep))/dTi), 1)*dTi];
istep = [istep; Ti-sum(istep)];

Tm  = 2000*year;
dTm = 20*year;
mstep = linspace(0.5*year, dTm, 5)';
mstep = [mstep; ones(floor((Tm-sum(mstep))/dTm),1)*dTm];
mstep = [mstep; Tm-sum(mstep)];

legendtext = {'No dissolution (A=0)', 'No dissolution (A=2)', ...
    'Dissolution (A=0)', 'Dissolution (A=2)'};
linetype = {'b--', 'r--', 'b-', 'r-'};
results = {};%cell(4,1);
fluid_types={'sharp interface','linear cap','P-scaled table','P-K-scaled table'};
% left and right panel of figure 15
%for n=1:2, % flat or non flat topsurface
residual=true;
pvt_types={'with_phase_boundary','coolprops_linear','coolprops_table'};

a=tic();

for n=1:2  
    figure(),clf,
    k=1;
    for i = 1:4 %residual saturation or not    
        %for i = 1:1 %residual saturation or not        
        %% Create model
        if n==1
            aquifer = makeAquiferModel('A',0);
            xc = aquifer.Gt.cells.centroids(:,1)/1e3;
            ff = 1;
        else
            aquifer = makeAquiferModel('A',2);
            xx = xc(150:650);
            ff = exp(-((xx-xc(400))/(0.3)).^2);
            ff = ff/sum(ff);                        
            % print -depsc2 figs/ex1-fig3a.eps;
            set(get(gca,'Children'),'LineStyle', '--');
        end
        
        G  = aquifer.G;
        Gt = aquifer.Gt;
        %% Make fluid model
        for dissolution=[false,true]            
        %for dissolution=[true]                
            fluid = makeFluidModel(aquifer, 'residual', residual, ...
                'dissolution', dissolution, 'fluidType', fluid_types{i});%,'co2_type',pvt_types{pvt_i});
            
            %% Setup system
            s = setupSimCompVe(aquifer.Gt, aquifer.rock2D);
            if(dissolution)
                systemOG = initADISystemVE({'Oil', 'Gas','DisGas'}, aquifer.Gt, aquifer.rock2D,...
                    fluid,'simComponents',s,'VE',true);
            else
                systemOG = initADISystemVE({'Oil','Gas'}, aquifer.Gt, aquifer.rock2D, ...
                    fluid, 'simComponents', s, 'VE', true);
            end
            systemOG.nonlinear.linesearch    = false;
            systemOG.nonlinear.maxIterations = 10;
            systemOG.nonlinear.tol           = 1e-6;
            %systemOG.nonlinear.tol           = 1e-5;
            
            %% Create well schedule
            z  = G.cells.centroids(:,3);
            clear x0;
            W  = aquifer.W;
            W(2).val = fluid.rhoOS * Gt.cells.z(W(2).cells)*norm(gravity);
            x0.pressure = W(2).val +(z(:)-z(W(2).cells))*norm(gravity)*fluid.rhoOS;
            x0.s(:,1)   = ones(G.cells.num,1);
            x0.s(:,2)   = zeros(G.cells.num,1);
            x0.rs       = ones(G.cells.num,1)*0.0;
            x0.smax     = x0.s;
            x0.smin     = x0.s;
            x0.sGmax    = x0.s(:,2);
            if(dissolution)
                x0.rs=ones(G.cells.num,1)*0.0;
            end
            
            control = struct('W',[],'step',struct('val',[],'control',[]));
            control.W = {W, W(2)};
            control.step.val = [istep; mstep];
            control.step.control = [ones(size(istep));ones(size(mstep))*2];
            
            %% Run the schedule setup
            t2=tic;
            [wellSols, states] = runMrstADI(x0, Gt, systemOG, control, ...
                'force_step', false, 'dt_min', 0.5*year, 'report_all', false);
            t2=toc(t2);
            xc=Gt.cells.centroids(:,1)/1e3;
            
            %% Plot results
            state = states{end-70};
            sG = free_sg(state.s(:,2),state.smax(:,2), ...
                struct('res_gas',fluid.res_gas, 'res_water', fluid.res_water));
            hold on
            %plot(xc, filter2(ff,sG.*Gt.columns.dz), linetype{k}, 'LineWidth', 2);
            hold off
            drawnow;
            
            results{end+1}=struct('states',{states}, 'ff', ff);%#ok
            k = k+1;
        end
    end
    axis tight
    set(gca,'YDir','reverse','FontSize',16);
    legend(legendtext{:}, 4);
    if(do_print)
        if(n==1)
            print -depsc2 figs/ex1-fig4a.eps;
        else
            print -depsc2 figs/ex1-fig4b.eps;
        end
    end
end
% print -depsc2 figs/ex1-fig3b.eps;
%save('data/disolutionTopSurfaceExample1Data','xc','results','control')
timing=toc(a);
save(['data/disolutionTopSurfaceExample1Data.mat'],'xc','results','control','timing');%#ok
%%
