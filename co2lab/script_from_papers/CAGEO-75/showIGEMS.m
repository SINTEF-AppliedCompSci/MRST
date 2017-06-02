%% Show the different types of surface morphologies for IGEMS
% The purpose of the IGEMS project was to investigate how variations in the
% top-surface morphology with a relief amplitude below seismic resolution
% influence CO2 storage capacity. To this end, the project considered a
% large 30x60 km sandbox in the shape of an inverted gutter. Fifteen
% different types of top-surface morphologies were designed by combining
% three different stratigraphic scenarios with five different structural
% scenarios. The stratigraphic scenarios are:
% 
% * flat deposition (flat)
% * burided beach ridges in a flooded marginal-marine setting (FMM),
% * buried offshore sand ridges (OSS).
% 
% The FMM deposition gives a dense set of small lobes with amplitude 1-10
% m, width 10-300 m, length less that 15 km, and spacing of 40-300 m. The
% OSS deposition gives significantly larger lobes with amplitude less than
% 20 meters, width 2-4 km, length 10-60 km, and spacing 2-4 km.
%
% The structural scenarios are:
%
% * no faults
% * uniform faults with a single 90 degree strike direction (UP1)
% * uniform faults with 30 and 90 degrees strike directions (UP2)
% * random faults with a single 90 degree strike direction (NP1)
% * random faults with 30 and 90 degrees strike directions (NP2)
%
% The fault displacement was 100 m in the uniform case and 20-150 m in the
% random case, and  the fault length was 4000 m in the uniform case and
% 300-6000 m in the random case.
%
% All data are publicly available and can be downloaded from the IGEMS
% project page <http://www.nr.no/IGEMS/>
%
% In this example, we compute fold and fault traps for one realization of
% each of the fifteen scenarios.

mrstModule add co2lab coarsegrid matlab_bgl;

%% Download data if necessary
idir = fullfile(mrstPath('co2lab'), 'data', 'igems');
if ~exist(fullfile(idir,'one_of_each'),'dir')
   disp(' -> Download data from: http://www.nr.no/IGEMS')
   disp(['    Putting data in ', idir]);
   unzip('https://www.nr.no/sites/default/files/files/one_of_each.zip', idir);
end


%% Display the data
idir = fullfile(idir,'one_of_each');
names = dir(idir);
names = names(3:end);
for i=1:numel(names)
   disp(['Realization: ' names(i).name]);
   Gt = topSurfaceGrid(...
      readIGEMSIRAP(names(i).name, [], 'dir', idir, 'save', false));
   
   res = trapAnalysis(Gt,true);
   z = 0*Gt.cells.z; 
   z(res.traps~=0) = res.trap_z(res.traps(res.traps~=0)) - Gt.cells.z(res.traps~=0);
   
   clf
   plotCellData(Gt, z, z>0,'EdgeColor','none');
   axis equal, axis([0 60000 0 30000]);
   view(-90,90);  box on, set(gca,'XTick',[],'YTick',[]); caxis([0 50]);
   % print('-dpng','-r300',[names(i).name(3:end-4) 'png']);
end