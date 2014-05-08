function param=run_opm_code_deck(deck)
%mrstModule add opm-utils
%mrstModule('add',fullfile(ROOTDIR,'mex','libgeometry'))
clear param;
param.use_deck='true';
param.output='true';
param.output_dir='test_output';
param.use_pside='false';
%param.linsolver='istl';
%param.linsolver='agmg';
param.linsolver='umfpack';
param.nl_pressure_masiter='10';
param.nl_pressure_change_tolerance='1';
param.nl_pressure_residual_tolerance='0';
param.nl_maxiter='30';
param.nl_tolerance='1e-4';
param.output_interval='1';
param.num_transport_substeps='1';
param.use_segregation_split='false';
writeDeck(deck,fullfile('data','test_deck'))
param.deck_filename=fullfile('data','test_deck','test_deck.DATA');
%% Run C OPM code
run_opm_core(param)
param=paramToStruct(fullfile(param.output_dir,'simulation.param'));

