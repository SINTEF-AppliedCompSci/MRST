%% make pie chart to show breakdown of trapping capacities of a formation
% then color the pie chart to show achieved storage of a particular well
% injection scenario and optimization.

% Inputs: fmCap.breakdowns, amounts achieved (computed from masses at last
% time step)


load coarsening_levels_70percent_of_full_StrapCap.mat;
n       = {names_and_cellsizes{:,1}};
c_level = {names_and_cellsizes{:,3}};


fmName = 'Tubaenfm';
[Gt, rock2D] = getFormationTopGrid(fmName, c_level{ find(strcmpi(fmName,n)) });

seainfo = getSeaInfo(fmName, 760);

% neglect dissolution trapping:
seainfo.dis_max = 0; % will remove an amount from strap and will remove btrap_diss

fmCap = getTrappingInfo(Gt, rock2D, seainfo, 'fmName', fmName);

%% Plot wedges of full storage potentials
labels = {['Struct. (',num2str(fmCap.breakdown.structural_trapping_capacity),' Gt)'],...
          ['Resid. (',num2str(fmCap.breakdown.residual_trapping_capacity),' Gt)'], ...
          ['Diss. (',num2str(fmCap.breakdown.dissolved_trapping_capacity),' Gt)']};
p = pie([fmCap.breakdown.structural_trapping_capacity, ...
    fmCap.breakdown.residual_trapping_capacity, ...
    fmCap.breakdown.dissolved_trapping_capacity], labels)

% modify colors of pie wedges
patches = findobj(p,'Type','patch');
patches(1).FaceColor = 'y'; % struct potential
patches(2).FaceColor = 'g'; % resid potential
patches(3).FaceColor = 'b'; % diss potential




%% Plot wedges of achieved masses and not acheived masses
totPotential = sum([fmCap.breakdown.structural_trapping_capacity, ...
    fmCap.breakdown.residual_trapping_capacity, ...
    fmCap.breakdown.dissolved_trapping_capacity]);
strapPerc = fmCap.breakdown.structural_trapping_capacity/totPotential;
rtrapPerc = fmCap.breakdown.residual_trapping_capacity/totPotential;
dtrapPerc = fmCap.breakdown.dissolved_trapping_capacity/totPotential;

strapAchieved = 1.62*0.412/totPotential;
rtrapAchieved = 4.13*0.029/totPotential;
dtrapAchieved = 0*0/totPotential;


figure;
p = pie([strapPerc - strapAchieved, strapAchieved, ...
    rtrapPerc - rtrapAchieved, rtrapAchieved, ...
    dtrapPerc - dtrapAchieved, dtrapAchieved], ...
    {'Structural','Achieved','Residual','Achieved','Dissolution','Achieved'});

% modify colors of pie wedges
patches = findobj(p,'Type','patch');

patches(1).FaceColor = 'y'; % struct not achieved
patches(1).FaceAlpha = 0.2;
patches(2).FaceColor = 'y'; % struct achieved

patches(3).FaceColor = 'g'; % resid not achieved
patches(3).FaceAlpha = 0.2;
patches(4).FaceColor = 'g'; % resid achieved

patches(5).FaceColor = 'b'; % diss not achieved
patches(5).FaceAlpha = 0.2;
patches(6).FaceColor = 'b'; % diss achieved

% shift pie wedges
for i = 1:2
    patches(i).XData = patches(i).XData - 0.1;
    patches(i).YData = patches(i).YData + 0.1;
end

for i = 3:4
    patches(i).XData = patches(i).XData - 0;
    patches(i).YData = patches(i).YData - 0.1;
end

for i = 5:6
    patches(i).XData = patches(i).XData + 0.1;
    patches(i).YData = patches(i).YData + 0.1;
end






