

%%% Compute volume of each trap in surface
names = [getBarentsSeaNames() getNorwegianSeaNames() getNorthSeaNames()];
% Remove certain formation names:
names = names(~strcmpi(names,'Nordmelafm'));
names = names(~strcmpi(names,'Rorfm'));
names = names(~strcmpi(names,'Notfm'));
names = names(~strcmpi(names,'Knurrfm'));       % @@ can be constructed??
names = names(~strcmpi(names,'Fruholmenfm'));   % @@
names = names(~strcmpi(names,'Cookfm'));
names = names(~strcmpi(names,'Dunlingp'));
names = names(~strcmpi(names,'Paleocene'));

N = 1;
figDirName = ['trap_vol_vs_height/ref' num2str(N)];
mkdir(figDirName);


for i=1:numel(names)
    
[Gt, rock2D] = getFormationTopGrid(names{i},N);
ta = trapAnalysis(Gt,false);

% %****
% tot_num_traps = max(unique(ta.traps(ta.traps~=0)));
% trap_vol = zeros(tot_num_traps,1);
% figure
% plotGrid(Gt, 'FaceColor','none', 'EdgeAlpha',0.1)
% for i=1:tot_num_traps
%     
%     plotCellData(Gt, Gt.cells.z, find(ta.traps == i))
%     
%     h = ta.trap_z(i) - Gt.cells.z(ta.traps == i); % @@ ?debug for cases where h is neg
%     v = h ...
%         .* Gt.cells.volumes(ta.traps == i);
%     trap_vol(i) = sum(v); % m3
%     clear v h
%     
% end
% %***see volumesOfTraps.m*** which gives different result than the above
% %code!

[trap_vol, h_max] = volumesOfTraps(Gt, ta); % total rock and poro vol, not pore volume.
% could use computeTrapVolume to get pore volume instead.

figure; set(gcf, 'name', names{i})
scatter(h_max, trap_vol, 'x', 'LineWidth',2)
xlabel('trap heights [m]','FontSize',16)
ylabel('trap volumes [m^3]','FontSize',16)
hax = gca; set(hax,'Fontsize',16);
box
drawnow

%saveas(gcf,[figDirName '/' names{i}],'fig')
export_fig(gcf, [figDirName '/' names{i}], '-png','-transparent')

clearvars -except i names figDirName N
close

end
