mrstModule add ad-core ad-blackoil deckformat diagnostics mrst-gui ad-props incomp optimization

 
%% Run EGG field simulation
Settup_Egg_simulation 
 


wellSols_ref =  wellSols;
model_ref    = model;
states_ref   = states;
schedule_ref = schedule;


%%  compute diagnostics 

 DD = WellPairNetwork(model_ref,schedule_ref,states_ref,state,wellSols_ref);
