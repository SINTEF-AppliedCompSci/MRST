% function d = setupEggRealizations

mrstModule add ad-blackoil deckformat ad-core diagnostics kpn-ds

% realizations = [2 13 5 17 22 2 56 5 1];
realizations = [1:100];

[models,wells] = getModelsEGG(realizations);
GlobalTT = cell(32,1);
GlobalPV = cell(32,1);
GlobalWellIndx = cell(32,1);

% simple diagnostics can be computed by calling
for  i =  1 : length(realizations)    
 [state, diagnostics] = computePressureAndDiagnostics(models{realizations(i)}, 'wells', wells{realizations(i)});
 
 %Calculate Well pair indices
 [IP_indices]=find(diagnostics.wellCommunication > 0);
 [I P]=find(diagnostics.wellCommunication > 0);
 %Calculate Transmisibility
     P_indx = P +8; %Producer index in wellsols

     for wp = 1:length(I) 
        fluxes(wp,1) = diagnostics.wellCommunication(I(wp),P(wp));

        DP(wp,1) = state.wellSol(I(wp)).pressure-... % Injector pressure
                   state.wellSol(P_indx(wp)).pressure;          %
        T(wp,1)  =   fluxes(wp,1)/DP(wp,1);
        pv(wp,1) = diagnostics.WP.vols(IP_indices(wp));
        
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
 

 
end

%% Post Procesing

for i =  1:32
    Mean_Trans(i) =  mean(GlobalTT{i});
    Min_Trans(i) =  min(GlobalTT{i});
    Max_Trans(i) =  max(GlobalTT{i});
    var_Trans(i) =  var(GlobalTT{i});
    count_Trans(i) =  length(GlobalTT{i});

    Mean_PV(i)   =  mean(GlobalPV{i});
    Min_PV(i)    =  min(GlobalPV{i});
    Max_PV(i)    =  max(GlobalPV{i});
    Var_PV(i)    =  var(GlobalPV{i});
    
    
    Mean_TT_PV(i)   =  mean(GlobalTT{i}./GlobalPV{i});
    Min_TT_PV(i)    =  min(GlobalTT{i}./GlobalPV{i});
    Max_TT_PV(i)    =  max(GlobalTT{i}./GlobalPV{i});
    Var_TT_PV(i)    =  var(GlobalTT{i}./GlobalPV{i});
    
    Well_Indices(i,:) = GlobalWellIndx{i}(1,:);
end

n_wells =numel(wells{1});
            for i = 1: n_wells
                cell_number = wells{1}(i).cells(1);
                XData(i) = models{1}.G.cells.centroids(cell_number,1);
                YData(i) = models{1}.G.cells.centroids(cell_number,2);
                ZData(i) = models{1}.G.cells.centroids(cell_number,3);
            end
            
   Graph = graph(Well_Indices(:,1),Well_Indices(:,2));         
            
figure            
subplot(2,3,1);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Mean_Trans/max(Mean_Trans))                        
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                        %plot(Graph,'k', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Min_Trans/max(Mean_Trans))
                        
                        %plot(Graph,'r','XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Max_Trans/max(Mean_Trans))

                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                title("Mean Transmisibility")
                
subplot(2,3,4);
                
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*var_Trans/max(var_Trans))                        
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                        %plot(Graph,'k', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Min_Trans/max(Mean_Trans))
                        
                        %plot(Graph,'r','XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Max_Trans/max(Mean_Trans))

                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                title("Variance Transmisibility")


 subplot(2,3,2);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'r', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Mean_PV/max(Mean_PV))
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ;                 
                 title("Mean Pore volume")            

subplot(2,3,5);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'r', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Var_PV/max(Var_PV))
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ;                 
                 title("Variance Pore volume")  
                 

 subplot(2,3,3);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'m', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Mean_TT_PV/max(Mean_TT_PV))
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ;                 
                 title("Mean (Transmisibility / Pore volume)")            

 subplot(2,3,6);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'m', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*Var_TT_PV/max(Var_TT_PV))
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ;                 
                 title("Variance (Transmisibility / Pore volume)") 
                 
                 
figure                                  
                subplot(1,2,1);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'k', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*count_Trans/100)
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                

                 title("More Frequent")         
 subplot(1,2,2);
                h2= plotGrid(models{1}.G, 'FaceColor', 'none', 'EdgeAlpha', 0.1), view(2);
                
                hold on, pg =  plot(Graph,'k', 'XData',XData,'YData',YData,'ZData',ZData,'LineWidth',10*(100-count_Trans+1)/100)
                        labelnode(pg,[1:n_wells],{wells{1}.name})
                hold off;
                % To have text indicating values with text                
                % plot(Graph,'EdgeLabel',compose("%5.2e",Graph.Edges.Weight),'LineWidth',10*TT/max(TT))
               
                pg.NodeFontSize= 15;
                axis off ; 
                

                 title("Less Frequent")            

                 
                

%[d] = DiagnosticsViewer(models,wells);