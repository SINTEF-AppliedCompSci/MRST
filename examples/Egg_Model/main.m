mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization

 
%% Run EGG field simulation
Settup_Egg_simulation 
 


wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;
W_ref        = schedule_ref.control.W;

%%  compute diagnostics 

 DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state,wellSols_ref);
 
 DD =  DD.filter_wps(5*stb/day);
 
 DD.plotWellPairConnections()
 
 
 
 % Initializing parameters
 for i =  1:numel(DD.wps)
    pv(i) = DD.wps{i}.volume;
    TT(i) = DD.wps{i}.Tr;
end

Tr_boxlimits = [0.1*TT ; 1.4*TT]';
Pv_boxlimits = [0.01*pv/10 ; 2*pv/10]';




%% Creatting data driven model

for i = 1:numel(W_ref)
wellname{i} = wellSols_ref{1}(1,i).name;
type{i} = wellSols_ref{1}(1,i).type;
end 
 wells.name = wellname;
 wells.type = type;           
 wells.val  =  mat2cell(vertcat(W_ref.val),ones(1,numel(W_ref)));
 
 
 
for i =  1:  numel(DD.wps)
    edges(i,:)= [DD.wps{i}.WellSolsIx_inj , DD.wps{i}.WellSolsIx_prod];
end            
    