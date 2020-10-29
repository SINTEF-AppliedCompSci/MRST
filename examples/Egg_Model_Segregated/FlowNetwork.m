mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization

 
%% Run EGG field simulation
Settup_Egg_simulation 
 


wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;

%%  Compute diagnostics 

 DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state,wellSols_ref);
 DD =  DD.filter_wps(0.05*stb/day);
 DD.plotWellPairConnections()
 DD.plotWellPairsData('subplot',[4,4])
 
 
 %% Making the Graph
 Well_Indices = DD.Graph.Edges.EndNodes;
 
 
 n_wells =numel(W_ref);
k = 1;
            for i = 1: n_wells
                Nodes(k,1) = k;
                Nodes(k+1,1) = k+1;
                Well(k,1) = i;
                Well_name{k,1}= [W_ref(i).name,'-',num2str(1)];
                Well(k+1,1) = i;
                cell_number = W_ref(i).cells(1);
                XData(k,1) = model.G.cells.centroids(cell_number,1);
                YData(k,1) = model.G.cells.centroids(cell_number,2);
                ZData(k,1) = model.G.cells.centroids(cell_number,3);
                Well_cell(k,1) = cell_number;
                SubWell(k,1) = 1;
                
                cell_number = W_ref(i).cells(end);
                Well_name{k+1,1}= [W_ref(i).name,'-',num2str(2)];
                XData(k+1,1) = model.G.cells.centroids(cell_number,1);
                YData(k+1,1) = model.G.cells.centroids(cell_number,2);
                ZData(k+1,1) = model.G.cells.centroids(cell_number,3);
                Well_cell(k+1,1) = cell_number;
                SubWell(k+1,1) = 2;
                k = k +2;
            end
            
            
            
            NodeWell_Indices = [];
            NodeWell_Indices_2 = [];
            k= 1;
            for i = 1: size(Well_Indices,1)
                    pv(k) = DD.wps{i}.volume;
                    TT(k) = DD.wps{i}.Tr;
                    pv(k+1) = DD.wps{i}.volume;
                    TT(k+1) = DD.wps{i}.Tr;
                    k=k+2;
                I1_well = find( Well_Indices(i,1) == Well);
                NodeWell_Indices = [NodeWell_Indices; I1_well];
                
                I2_well = find( Well_Indices(i,2) == Well);
                NodeWell_Indices_2 = [NodeWell_Indices_2; I2_well];
            end
            
        EdgeTable = table([NodeWell_Indices NodeWell_Indices_2],'VariableNames',{'EndNodes'});
        NodeTable = table(Nodes,Well,SubWell,Well_name,Well_cell,XData,YData,ZData,'VariableNames',{'Nodes','Well','SubWell','Well_name','Well_cell','XData','YData','ZData'});
        Graph = graph(EdgeTable,NodeTable);            

   

%% Creatting data driven model

L = 435;
G = cartGrid([10, 1, numedges(Graph)], [L, L/5 ,L/5]*meter^3);
G = computeGeometry(G);


fluid =  model_ref.fluid;
rock = makeRock(G, 1000*milli*darcy, 0.1);

gravity off
model = GenericBlackOilModel(G, rock, fluid);
model.gas=false;
model.OutputStateFunctions = {};

 
%[model,W,indexs] = createDDmodel_1(model,10,DD.Graph,W_ref);

obj = GPSNet(model,10,Graph,W_ref);
model = obj.model;
W     = obj.W;
indexs.faces = obj.Graph.Edges.Face_Indices;
indexs.cells = obj.Graph.Edges.Cell_Indices;

figure
subplot(1,2,1)
obj.plotNetwork(model_ref.G,pv)

subplot(1,2,2)
obj.plotNetwork(model_ref.G,TT)


figure
obj.plotGPSNetGrid([])