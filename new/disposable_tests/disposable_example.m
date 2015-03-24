function disposable_example()
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

aquifer = makeAquiferModel('A',0,'D',depth);
xc      = aquifer.Gt.cells.centroids(:,1)/1e3;
ff      = 1;
            
Gt = aquifer.Gt;

fluid = makeVEFluid(Gt, aquifer.rock2D, 'simple');

model = CO2VEBlackOilTypeModel(Gt, aquifer.rock2D, fluid);





end
         