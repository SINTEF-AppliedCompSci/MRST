function deckmrst = modelToDeckSections(model, state0, schedulemrst)
nc=model.G.cells.num;
nif=size(model.operators.N(:,1),1);

RUNSPEC=[];%deck.RUNSPEC;
RUNSPEC.DIMENS=[model.G.cells.num,1,1];

%SATOPTS=[]
%SATOPTS.DIRECT=0;
%SATOPTS.HYSTER=0;
%SATOPTS.IRREVERS=0;
%RUNSPEC.SATOPTS=SATOPTS;
%RUNSPEC.GRIDOPTS={'YES'  [0]  [0]};
RUNSPEC.OIL=model.oil;
RUNSPEC.GAS=model.gas;
RUNSPEC.WATER=model.water;
RUNSPEC.DISGAS=model.disgas;
RUNSPEC.VAPOIL=model.vapoil;
RUNSPEC.UNIFOUT= 1;
RUNSPEC.METRIC=1;
RUNSPEC.WELLDIMS = [10 3 10 20 5 10 5 4 3 0 1 1];
%RUNSPEC.REGDIMS= [1 1 1 0 0 1 0 0 0];
%RUNSPEC.EQLDIMS= [1 100 50 1 50];
RUNSPEC.START=734813;
RUNSPEC.TABDIMS=[1 1 33 60 16 60 20 1 1 1 10 1 -1 0 1];
% model to grid

GRID=[];
GRID.INIT=1;
%dx=(model.G.cells.volumes).^(1/3);

GRID.DXV=model.G.cells.volumes;
GRID.DYV=1;
GRID.DZV=2*1e3;
%GRID.TOPS=dx*0;
%GRID.PERMX=ones(nc,1);
%GRID.PERMY=ones(nc,1);
%GRID.PERMZ=ones(nc,1);
%GRID.PORO
GRID.ACTNUM=int32(ones(nc,1));
GRID.cartDims=RUNSPEC.DIMENS;
T = model.operators.T;
T = T/((centi*poise*meter^3)/(day*barsa));
GRID.NNC=[model.operators.N(:,1),ones(nif,2),model.operators.N(:,2),ones(nif,2), T];
% model to edit
EDIT.PORV=model.operators.pv;
EDIT.DEPTH=model.G.cells.centroids(:,3);



%%
REGIONS=struct();
SUMMARY=struct();
% SOLUTION
SOLUTION.PRESSURE=state0.pressure/barsa;
SOLUTION.SWAT=state0.s(:,1);
SOLUTION.SGAS=state0.s(:,3);
SOLUTION.SOIL=state0.s(:,2);
SOLUTION.RS=state0.rs;
SOLUTION.RV=state0.rv;
%%
UnhandledKeywords=[];
UnhandledKeywords.SUMMARY = {'ALL'};
UnhandledKeywords.SCHEDULE= {'OPTIONS'  'RPTRST'};

%deckmrst=deck;
deckmrst=[];
deckmrst.GRID=GRID;
deckmrst.RUNSPEC=RUNSPEC;
deckmrst.EDIT=EDIT;
deckmrst.REGIONS=REGIONS;
deckmrst.SUMMARY=SUMMARY;
deckmrst.SOLUTION=SOLUTION;
deckmrst.UnhandledKeywords=UnhandledKeywords;
%%
%deckmrst_unstr.PROPS=deck_mrst_org.PROPS;% need to be implmented
%deckmrst_unstr.SCHEDULE=deck_mrst_org.SCHEDULE;
deckmrst.RUNSPEC.UNIFOUT=1;
deckmrst.GRID.INIT=1;
deckmrst.SCHEDULE.RPTSCHED=1;
G_unstr=model.G;
G_unstr.type={'unstr'};
G_unstr.cartDims=[G_unstr.cells.num,1,1];
G_unstr.cells.indexMap=1:G_unstr.cells.num;
for i=1:numel(schedulemrst.control)
    deckmrst.SCHEDULE.control(i)=mrstWellToControl(schedulemrst.control(i).W, G_unstr, struct() ,deckmrst.RUNSPEC,'add_wellindex',true);
end
deckmrst.SCHEDULE.step = schedulemrst.step;
deckmrst.SCHEDULE.step.val = deckmrst.SCHEDULE.step.val/day;