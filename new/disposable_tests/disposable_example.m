function [Gt, wellSols, states] =  disposable_example()

   gravity on;
   moduleCheck('ad-core');
% Example constructed to test the new framework.  Throw away when things work

Ti      =   50*year;
dTi     =  2*year;
istep   = linspace(0.1*year, dTi, 10)';
istep   = [istep; ones(floor((Ti-sum(istep))/dTi), 1)*dTi];
istep   = [istep; Ti-sum(istep)];

Tm      = 2000*year;
dTm     = 20*year;
mstep   = linspace(0.5*year, dTm, 5)';
mstep   = [mstep; ones(floor((Tm-sum(mstep))/dTm),1)*dTm];
mstep   = [mstep; Tm-sum(mstep)];
depth   = 1300;

aquifer = makeAquiferModel('A',0,'D',depth);
xc      = aquifer.Gt.cells.centroids(:,1)/1e3;
ff      = 1;
            
% Adapting wells to new structure
aquifer.W(1).compi = [0 1];
aquifer.W(2).compi = [1 0]; % we now treat brine as the _water_ component, not oil


Gt = aquifer.Gt;

fluid = makeVEFluid(Gt, aquifer.rock2D, 'simple');

model = CO2VEBlackOilTypeModel(Gt, aquifer.rock2D, fluid);

% Setup state
bhp = aquifer.W(2).val;
z   = aquifer.G.cells.centroids(:,3);
initState.pressure = bhp + (z(:) - z(aquifer.W(2).cells)) * norm(gravity) * fluid.rhoWS;
initState.s = [ones(Gt.cells.num, 1), zeros(Gt.cells.num, 1)];
initState.sGmax = initState.s(:,2);

% Setup schedule
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

Wctrl1 = aquifer.W;
Wctrl2 = Wctrl1;
Wctrl2(1).val = 0;
schedule = struct('control', [struct('W', Wctrl1), struct('W', Wctrl2)], ...
                  'step', struct('control', [ones(size(istep)); 2 * ones(size(mstep))], ...
                                 'val', [istep; mstep]));
               
% Run simulation
[wellSols, states] = simulateScheduleAD(initState, model, schedule);

end
         