%%
%dir='/home/hnil/Documents/data/benchmark_new'
mrstModule add ad-core ad-blackoil
dir='/home/hnil/Documents/GITHUB/opm-data/norne'
ifile='NORNE_ATW2013.DATA'
deck_file=fullfile(dir,ifile);
mrstModule add deckformat
deck  = readEclipseDeck(deck_file);
deck = convertDeckUnits(deck);
G=initEclipseGrid(deck);
G=G(1);
%%
G=computeGeometry(G)
fluid=initEclipseFluid(deck);
rock=initEclipseRock(deck);
rock=compressRock(rock,G.cells.indexMap);
model = selectModelFromDeck(G, rock, fluid, deck);
%schedule = convertDeckScheduleToMRST(model, deck,'strictParsing',false);
W=processWells(model.G, model.rock, deck.SCHEDULE.control(end),'strictParsing',false)

%%'strictParsing',false
mrstModule add vemmech opm_gridprocessing%project-mechanics-fractures 
%%
% make extendend grid
grd=deck.GRID;
grd.ACTNUM=ones(size(grd.ACTNUM));
actnum=reshape(grd.ACTNUM,grd.cartDims);
grd.ACTNUM=reshape(actnum,[],1);
grdc = verticalGrdecl(grd);
% pad grid on top and bottom
grdc=padGrdecl(grdc,[0 0 1],[0 10000;0 5000;[1000 500]/50]);
grdc = verticalGrdecl(grdc);
ref=4; %refine padding
grdc = refineGrdeclLayers(grdc, [1 1], ref);
grdc = refineGrdeclLayers(grdc, [grdc.cartDims(3) grdc.cartDims(3)], ref);
G_cut=processGRDECL(grdc,'tolerance',1);%,'RepairZCORN',true,'Tolerance',1);

clf,
%plot padded grid
plotGrid(G_cut),view(3)
%%
%mrstModule add vemmech
%%
% calculate geometry and extra mappings needed for mechanics
profile off;profile on
G_cut=computeGeometry(G_cut);
Ga = createAugmentedGrid(G_cut);
profile off;profile viewer
%%
% set bounary conditions
zmax = max(Ga.nodes.coords(:,3));
nodes = find(zmax==Ga.nodes.coords(:,3));
mask = true(numel(nodes),Ga.griddim);
disp_node =zeros(numel(nodes),Ga.griddim);
disp_bc = struct('nodes', nodes, 'uu', disp_node, 'faces', [], 'uu_face', [], 'mask', mask);
force_bc=[];
%% 
%% read in simulation from eclipse on the norne case
ofilen= 'NORNE_ATW2013';
ecl_case=false;
if(ecl_case)
    odir = '/home/hnil/Documents/sims/opm/norne_ecl';
    %ofilen= 'NORNE_ATW2013';
    fname=fullfile(odir,ofilen);
    out=readEclipseRestartUnFmt(fname);
else
%%
    mrstModule add project-flowdiagnostics-resinsight
    mrstModule add ad-fi % to be able to use functions moved to private in deckformat
    mrstModule add deckformat % to avoid using old versions
    odir = '/home/hnil/Documents/sims/opm/norne_flow';
    fname=fullfile(odir,ofilen);
    out=readRestart(fname);
end
%
ifile=fullfile(odir,[ofilen,'.INIT']);
init=readEclipseOutputFileUnFmt(ifile);
if(~ecl_case)
   init.NTG.values = []; % for compatibilty 
end
egridf=fullfile(odir,[ofilen,'.EGRID']);
egrid=readEclipseOutputFileUnFmt(egridf);
% make grid from the output files of the norne case 
[Ge,rock,T] = eclOut2mrst(init,egrid);
%%
myfn =@(fig) set(fig,'FontSize',24)
myprint =@(name) print('-dpng',['figs/norne_',name,'.png']);
%myprint =@(name) disp('noprint');
myview =@() view([1 -4 10]); 
mypos  =@(fig) set(fig,'Position',[200 230 1050 850]);
opts_box={'EdgeAlpha',0.2};
figure(11),clf
plotCellData(Ge,out.PRESSURE{end}-out.PRESSURE{1}),colorbar,axis off
plotWell(G,W)
myview();
mypos(gcf);
myfn(gca)
myprint('pressure_diff')
ca=caxis();
figure(21),clf
colorbar
caxis(ca),axis off,
myfn(gca)
myprint('pressure_diff_col')
%%
figure(12),clf
plotCellData(Ge,out.SWAT{end}),colorbar,axis off
plotWell(G,W)
myview();
mypos(gcf)
myprint('SWAT_end')
%%
figure(13),clf
plotCellData(Ge,out.SWAT{1}),colorbar,axis off
plotWell(G,W)
myview();
mypos(gcf)
myprint('SWAT_1')
%%
figure(14),clf
plotCellData(Ge,out.PRESSURE{end}),colorbar,axis off
plotWell(G,W)
myview();
mypos(gcf)
myprint('pressure_end')

%%
% set forces for linear ealsticity to differences in pressure
shift=prod(Ga.cartDims(1:2))*ref;
press=zeros(prod(Ga.cartDims),1);
press(Ge.cells.indexMap+shift)=out.PRESSURE{end}-out.PRESSURE{1};
pressure=press(Ga.cells.indexMap);
figure(33),clf
plotGrid(Ga,'FaceColor','none',opts_box{:})
plotCellData(Ga,pressure,find(abs(pressure)>0)),colorbar,axis off,
myview()
mypos(gcf)
myprint('press_diff_extended')
%%
figure(44),clf
bc_fixed=pside([],Ga,'BOTTOM',0)
bc_o=pside([],Ga,'LEFT',0)
bc_o=pside(bc_o,Ga,'RIGHT',0)
bc_o=pside(bc_o,Ga,'TOP',0)

plotGrid(Ga,'FaceColor','none',opts_box{:})
plotCellData(Ga,pressure,find(abs(pressure)>0))%,colorbar
axis off,
plotFaces(Ga,bc_fixed.face,'FaceColor','g')
plotFaces(Ga,bc_o.face,'FaceColor','m','FaceAlpha',0.1)
myview()
mypos(gcf)
myprint('press_diff_extended_bc')
%%
% put material popropertis
opt=struct('E',1e9,'nu',0.3);
Ev     = repmat(opt.E, Ga.cells.num, 1);
nuv    = repmat(opt.nu, Ga.cells.num, 1);
C      = Enu2C(Ev, nuv, Ga); % full hooks tensor

%final boundary structure
el_bc = struct('disp_bc', disp_bc, 'force_bc', force_bc);
lload=@(x) 0*x;
bbsize = 30000-(Ga.griddim-2)*20000;
% SOLVE mechanics
lsolve = @mldivide;
profile off;profile on;
[uu,extra] = VEM_linElast(Ga, C, el_bc, lload, 'linsolve', lsolve, 'blocksize', bbsize, ...
                  'force_method', 'dual_grad_type','pressure',pressure*barsa,'experimental_scaling', false) ;
profile off;profile viewer
%save('none_mech','uu','extra','Ga','-v7.3')
%save('none_mech_uu','uu','-v7.3')
%save('none_mech_uu_new','uu','-v7.3')
%{
[uu,extra] = VEM_linElast(Ga, C, el_bc, lload, 'linsolve', lsolve, 'blocksize', bbsize, ...
                  'force_method', 'dual_grad_type','pressure',pressure*barsa,'experimental_scaling', false,'no_solve',true) ;
a=load('none_mech_uu_new');
uu=a.uu;
%}
%% Pot displace ment on embeded grid 
cells=find(abs(pressure)>0); %cells of embedded grid
figure(15)
clf,
plotGrid(Ga,'FaceColor','none',opts_box{:})
% plot displace ment in z direction
plotNodeData(Ga,uu(:,3),'cells',cells),colorbar,
myview(),axis off
mypos(gcf)
myprint('dz_embedded');
%%
%% Pot displace ment on embeded grid 
cells=find(abs(pressure)>0); %cells of embedded grid
figure(16)
clf,
%plotGrid(Ga,'FaceColor','none','EdgeAlpha',0.3)
% plot displace ment in z direction
plotNodeData(Ga,uu(:,3)),colorbar,
myview(),axis off
mypos(gcf)
myprint('dz_outer');
%%
div=VEM3D_div(Ga);
%%
myfont =@(fig) set(fig,'FontSize',16)
figure(17),clf
%plotGrid(Ga,'FaceColor','none',opts_box{:})
plotCellData(Ga,div*reshape(uu',[],1)./Ga.cells.volumes,cells);colorbar
myview(),axis off
mypos(gcf)
myfont(gca)
myprint('div_norne');
%%
% calculate stress changes
uue=uu-1e-1*[zeros(Ga.nodes.num,2) Ga.nodes.coords(:,3)];
stress=calculateStressVEM(Ga,uue,extra);
%calculate eigen values of stress changes
[sigm,evec]=calStressEigsVEM(Ga,stress);
% plot smalles eigenvalue of stress chances
%%
figure(17)
clf,
plotCellData(Ga,sigm(:,1),cells);colorbar,axis off
myview()
mypos(gcf)
myprint('simg1_norne');
%%
figure(18)
clf,
plotCellData(Ga,sqrt(evec(:,1).^2+evec(:,2).^2),cells);colorbar,axis off
%plotCellData(Ga,sqrt(evec(:,1).^2+evec(:,2).^2),cells);colorbar,axis off
myview()
mypos(gcf)
myprint('vec12_norne');

%caxis([min(sigm(cells,1)),3.1e7])
%caxis([4e7 4.5e7])
%clf,plotCellData(Ga,evec(:,2),cells);colorbar
%%
% here one can test solvers

%%
A=extra.A;b=extra.rhs;
%profile off;profile on
%%
disp('mldivide')
tic;
%linsolve=@mldivide;
%linsolve=@pcg(A,rhs,
x=A\b;
toc;
%%
%%{
%disp('pcg ichol')
mrstModule add project-mechanics-fractures
profile off;profile on
At=(A+A')/2;
bt=b;
L = shiftedIChol(At, 0.5, 'droptol', 1e-1, 'type', 'nofill');
x1 = itersolve_JP(At, bt, 'solver', @pcg, ...
                             'M1', L, 'M2', L', 'x0', zeros(size(x)), 'tol', 1e-4);
profile off;profile viewer
%toc;
%%                         
%L = ichol(A,struct('type','ict','droptol',1e-03,'michol','on'));
%opts.type = 'nofill';
%opts.michol = 'on';
%tol=1e-5;maxit=10;
%L = ichol(A,opts);
%[x1,flag1,rr1,iter1,rv1] = pcg(A,b,tol,maxit,L,L');
toc;
%}
%%
%%{
disp('bicg ilu')
tic;
setup.type = 'nofill';
setup.milu = 'row';
setup.droptol = 0.1;
[L,U] = ilu(A,setup);
toc;
[x1,flag1,rr1,iter1,rv1] = pcg(A,b,1e-5,40,L,U);
%}
%profile off;profile viewer

%%
% <html>
% <p><font size="-1">
% Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
