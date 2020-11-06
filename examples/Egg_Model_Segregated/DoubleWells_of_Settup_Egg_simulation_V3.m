%% Example demonstrating the two-phase oil-water Egg model
% This example sets up and runs the Egg model using the two-phase AD
% solvers. 
%
% For details on the EggModel and the corresponding ensamble, see
% Jansen, J. D., et al. "The egg model–a geological ensemble for reservoir
% simulation." Geoscience Data Journal 1.2 (2014): 192-195.

mrstModule add ad-core ad-blackoil deckformat diagnostics

% Realizations can be set to 0 for base cae, or a number between 1 and 100
% for different permeabilities.sa
%realization = [2];
[G, rock, fluid, deck] = setupEGG();

%% Changing the well cells for producer form 1:7 to 1:3
for i = 9:12
    deck.SCHEDULE.control.COMPDAT{i,5} = [3];
end

[state, model, schedule, nonlinear] = initEclipseProblemAD(deck, 'G', G, 'TimestepStrategy', 'none');

load('Egg_Rock.mat','rock')

model.getPhaseNames()

state=initResSol(model.G, 200*barsa);
hT = computeTrans(model.G, model.rock);
state = incompTPFA(state, model.G, hT, model.fluid,'wells', schedule.control.W);
%% Adding oil top and water bottom
for i = 1:model.G.cells.num
    if model.G.cells.centroids(i,3) < 4014
        state.s(i,1:2) = [0 1];
    else % model.G.cells.centroids(i,3) >= 4010
        state.s(i,1:2) = [1 0];
    end
end

% plotGrid(G, 'FaceColor', 'none', 'EdgeAlpha', 0.1)
% plotWell(G,schedule.control.W)


problem = packSimulationProblem(state, model, schedule, 'EGG_realization_0_segregated', 'NonLinearSolver', nonlinear);

%%
[ok, status] = simulatePackedProblem(problem);


%% Run simulation
[wellSols, states, reports] = getPackedSimulatorOutput(problem);
    
  %d = PostProcessDiagnosticsMRST(problem);
 
 %state= states{120};
 W = schedule.control.W;
 nw = numel(W);
 nsegment = 2;
 [statexd,Wdf] = expandWellCompletions(state, W,[(1:nw)' repmat(nsegment,nw,1)]);
 

 
 
 %% Flow diagnostics
    Df  = computeTOFandTracer(statexd, model.G, model.rock, 'wells', Wdf);
    WPf = computeWellPairs(statexd, model.G, model.rock, Wdf, Df);
 salloc = cellfun(@sum, {WPf.inj.alloc}, 'UniformOutput',false);
 wellCommunication = vertcat(salloc{:});

 [d] = DiagnosticsViewer({model},{Wdf},'state0',{statexd});
  %%  PostProcesinf Initial Saturation

 D = d.Data{1}.diagnostics.D;
 
 % For each wellpair
 i_producer = 1;  % P1.1 
 
    btof = D.ptof(:,i_producer);
    s = d.Data{1}.states.s(:,1);
    [btof, ix] = sort(btof);
    s = s(ix);
    
    pv = d.Data{1}.static(6).values;
    pv=pv(ix);
    target =cumsum(s.*pv)./cumsum(pv);
    
    NN= 6610;
    indexes = 1:100:NN;
    
    x = pv(indexes);
    y = target(indexes);%s(indexes);

    yy = spline(x,y,btof(1:NN));
    
    
    
    figure    
    subplot(1,2,1);
    plot(btof/year,target,x/year,y,'or');
    legend('Cumulative Saturation', 'linear interpolation')
    
    %derivate = diff(yy)./diff(pv(1:NN));
    subplot(1,2,2);
    plot(btof(1:NN-1)/year,s(1:NN-1),'o');
    
 %%  PostProcesinf FlowDiagnostics

 nit   = numel(Df.inj);
 npt   = numel(Df.prod);
 nseg  = nit + npt;
  cmap = lines(nseg); %cmap = 0.6*cmap + .4*ones(size(cmap));
 colormap(cmap);
 figure; clf
    %set(gcf,'Position', [860 450 840 310],'PaperPositionMode','auto');
    subplot(1,2,1);
    for i=1:nit
       plotCellData(G,Df.ipart,Df.ipart==i,'FaceAlpha',.3,'EdgeAlpha',.2,'FaceColor',cmap(i,:));
    end
    plotGrid(G,'FaceColor','none');
    for i=1:nseg
       plotGrid(G,Wdf(i).cells,'FaceColor',cmap(i,:));
    end
    plotWell(G,W,'Color','k');
    %view(-40,10); 
    axis tight off; caxis([1 nit+npt]);
    
    subplot(1,2,2);

    for i=1:npt
      plotCellData(G,Df.ppart+nit,Df.ppart==i,'FaceAlpha',.4,'EdgeAlpha',.2,'FaceColor',cmap(i,:));
    end
    plotGrid(G,'FaceColor','none');
    for i=1:nseg
       plotGrid(G,Wdf(i).cells,'FaceColor',cmap(i,:));
    end
    plotWell(G,W,'Color','k');
    %view(-40,10); 
    axis tight off; caxis([1 nit+npt]);
    colormap(cmap);
 
 %% Procesing 
 %Calculate Well pair indices
 [IP_indices]=find(wellCommunication > 0);
 [I P]=find(wellCommunication > 0);
 %Calculate Transmisibility
     P_indx = P +16; %Producer index in wellsols
GlobalTT = cell(128,1);
GlobalPV = cell(128,1);
GlobalWellIndx = cell(128,1);

     for wp = 1:length(I) 
         wp
        fluxes(wp,1) = wellCommunication(I(wp),P(wp));

        DP(wp,1) = statexd.wellSol(I(wp)).bhp-... % Injector pressure
                   statexd.wellSol(P_indx(wp)).bhp;          %
        T(wp,1)  =   fluxes(wp,1)/DP(wp,1);
        pv(wp,1) = WPf.vols(IP_indices(wp));
        
        GlobalTT{IP_indices(wp)} =  [GlobalTT{IP_indices(wp)},T(wp,1)];
        GlobalPV{IP_indices(wp)} =  [GlobalPV{IP_indices(wp)},pv(wp,1)];
        GlobalWellIndx{IP_indices(wp)} =  [GlobalWellIndx{IP_indices(wp)};[I(wp),P_indx(wp)]];

     end
     
% 
%      EndNodes = [I,P_indx];
%      Counter  = 0*T;
%      EdgeTable = table([I P_indx],T,pv,fluxes,DP,Counter,'VariableNames',{'EndNodes' 'T' 'pv' 'fluxes' 'DP' 'Counter'});
%      %Indicate a counter on that edge
%      
%  Graph = graph(EdgeTable);
%% Post Procesing
k = 1;
for i =  1:128
    if ~isempty(GlobalTT{i})
        
        Mean_Trans(k) =  mean(GlobalTT{i});
        Min_Trans(k) =  min(GlobalTT{i});
        Max_Trans(k) =  max(GlobalTT{i});
        var_Trans(k) =  var(GlobalTT{i});
        count_Trans(k) =  length(GlobalTT{i});

        Mean_PV(k)   =  mean(GlobalPV{i});
        Min_PV(k)    =  min(GlobalPV{i});
        Max_PV(k)    =  max(GlobalPV{i});
        Var_PV(k)    =  var(GlobalPV{i});


        Mean_TT_PV(k)   =  mean(GlobalTT{i}./GlobalPV{i});
        Min_TT_PV(k)    =  min(GlobalTT{i}./GlobalPV{i});
        Max_TT_PV(k)    =  max(GlobalTT{i}./GlobalPV{i});
        Var_TT_PV(k)    =  var(GlobalTT{i}./GlobalPV{i});

        Well_Indices(k,:) = GlobalWellIndx{i}(1,:);
        k = k+1;
    end
end
%% Creating the Graph

% Nodes

n_wells =numel(Wdf);
            k =1;
            for i = 1:numel(W)
                Well(k,1) = i;
                Well(k+1,1) = i;
                SubWell(k,1) = 1;
                SubWell(k+1,1) = 2;
                k = k+2;
            end

            for i = 1: n_wells
                Nodes(i,1) =i;
                Well_name{i,1}= Wdf(i).name;
                cell_number = Wdf(i).cells(1);                
                XData(i,1) = model.G.cells.centroids(cell_number,1);
                YData(i,1) = model.G.cells.centroids(cell_number,2);
                ZData(i,1) = model.G.cells.centroids(cell_number,3);
                Well_cell(i,1) = cell_number;                
            end


        EdgeTable = table(Well_Indices,'VariableNames',{'EndNodes'});                   
        NodeTable = table(Nodes,Well,SubWell,Well_name,Well_cell,XData,YData,ZData,'VariableNames',{'Nodes','Well','SubWell','Well_name','Well_cell','XData','YData','ZData'});
        Graph = graph(EdgeTable,NodeTable); 
        
        
        