load ib
load id
load rb
load rd

data = { {ib, rb}, {id, rd} };
wname = { 'B-2H', 'D-3AH' };
fname = { 'Time [days]', 'Bottom hole pressure [bars]', 'Oil rate [m^3/day]', 'Gas rate [m^3/day]', 'Water rate [m^3/day]' };
sname = {'time', 'pressure', 'orat', 'grat', 'wrat'};


close all
fs = 20;

pth = fullfile(mrstPath('dg'), 'examples', 'ecmor-xvi', 'norne', 'fig');
saveeps = @(name) print(fullfile(pth, name), '-depsc');

for well = 1:2
    wdata = data{well};
    imp = wdata{1};
    reo = wdata{2};
    for field = 2:5
        figure
        hold on
        plot(imp(:,1), imp(:,field), '-', 'linew', 2)
        plot(reo(:,1), reo(:,field), 'o', 'linew', 2)
        hold off
        if well == 1 && field == 2
            legend('Fully implicit', 'Reordering');
        end
%         xlabel(fname{1})
%         ylabel(fname{field})
        xlim([0, max(imp(:,1))]);
        ax = gca;
        ax.FontSize = fs;
        box on
        drawnow()
        pause(0.1)
        saveeps([wname{well}, '-', sname{field}])
%         title([fname{field} ' of well ' wname{well}])
    end
end