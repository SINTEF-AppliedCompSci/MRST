mrstModule add opm-utils
mrstModule('add',fullfile(ROOTDIR,'mex','libgeometry'))
mrstModule add spe10
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
%layers=[1:85];
layers=[1];
deck=deckQFS(layers);
nstep=40;
deck.SCHEDULE.step.control=ones(nstep,1);
deck.SCHEDULE.step.val=ones(nstep,1)*2000/nstep;
deck.PROPS.PVCDO=deck.PROPS.PVTW;
deck.PROPS.PVCDO(4)=deck.PROPS.PVCDO(4)*10;
deck.PROPS=rmfield(deck.PROPS,'PVDO');
s=linspace(0,1,10)';alpha=2;
deck.PROPS.SWOF{1}=[s,s.^alpha,(1-s).^alpha,s*0];
deck.PROPS.ROCK=[0 0 NaN NaN NaN NaN]
deck.GRID.PORO(deck.GRID.PORO==0)=min(deck.GRID.PORO(deck.GRID.PORO>0));
initsat=0.0;
deck.SOLUTION.SOIL(:)=1-initsat;
deck.SOLUTION.SWAT(:)=initsat;
writeDeck(deck,fullfile('data','test_deck'))
param.deck_filename=fullfile('data','test_deck','test_deck.DATA');
%%

output=run_opm_core(param)
%% read deck used in simulation
%param=paramToStruct(fullfile(param.output_dir,'simulation.param'));
deck=readEclipseDeck(param.deck_filename);
deck_orig=deck;
deck.GRID.ACTNUM=int32(ones(prod(deck.GRID.cartDims),1));
deck=convertDeckUnits(deck);
%% init grid
init_model=true;
if(init_model)
   g=initEclipseGrid(deck);
   g=mcomputeGeometry(g);
   rock=initEclipseRock(deck);
   rock = compressRock(rock, g.cells.indexMap);
end
%% read results
nsteps=numberOfsteps(param.output_dir)-1;
opmres = getOpmResults(param.output_dir, nsteps);
%{
press=load(fullfile(param.output_dir,'pressure-001.dat'))/barsa;
sat=load(fullfile(param.output_dir,'saturation-001.dat'));
figure(1),clf
plotCellData(g,press);colorbar;caxis([min(press) max(press)])
figure(2),clf
sat=reshape(sat,2,[])';
plotCellData(g,sat(:,1),sat(:,1)>0.1);colorbar;caxis([0 1])
%}
press=opmres(end).pressure/barsa;
sat=opmres(end).saturation;
figure(4),clf
subplot(2,1,1),cla
plotCellData(g,press);colorbar;caxis([min(press) max(press)]);axis tight
subplot(2,1,2),cla
plotCellData(g,sat(:,1),sat(:,1)>0.1);colorbar;caxis([0 1]);axis tight
wells=load('test_output/wellreport.txt');
figure(3)
wells(:,4)=1;
plot(wells(:,1),wells(:,4:3:end))
for i=1:numel(W)
   [d,ind]=min(abs(wells(:,4+3*(i-1))-((wells(end,1)-wells(:,1))./wells(end,1))));
   text(wells(ind,1),wells(ind,4+3*(i-1))-0.05,W(i).name)
end
figure(1)
for i=1:3
   subplot(3,1,i)
   plot(wells(2:end,1),wells(2:end,1+i:3:end))
end
