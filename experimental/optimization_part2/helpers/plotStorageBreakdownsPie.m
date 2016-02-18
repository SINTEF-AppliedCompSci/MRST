function plotStorageBreakdownsPie( strap_pot, rtrap_pot, dtrap_pot, strap_ach, rtrap_ach, dtrap_ach )
% Makes pie chart of trapping breakdowns, both full potential and achieved
% Ensure units are consistently in Gt.

    % Adjust shifting of pie wedges according to slice sizes.
    % Deal with no dissolution trapping.

    %% Plot wedges of full storage potentials
    labels = {['Struct. (',num2str(strap_pot),' Gt)'],...
              ['Resid. (',num2str(rtrap_pot),' Gt)'], ...
              ['Diss. (',num2str(dtrap_pot),' Gt)']};
    p = pie([strap_pot, rtrap_pot, dtrap_pot], labels);

    % modify colors of pie wedges
    patches = findobj(p,'Type','patch');
    patches(1).FaceColor = 'y'; % struct potential
    patches(2).FaceColor = 'g'; % resid potential
    patches(3).FaceColor = 'b'; % diss potential

    %% Plot wedges of achieved masses and not acheived masses
    totPotential = sum([strap_pot, rtrap_pot, dtrap_pot]); % Gt
    strapPerc = strap_pot/totPotential * 100; % percent of total
    rtrapPerc = rtrap_pot/totPotential * 100;
    dtrapPerc = dtrap_pot/totPotential * 100;

    strapAchieved = strap_ach/totPotential * 100; % percent of total
    rtrapAchieved = rtrap_ach/totPotential * 100;
    dtrapAchieved = dtrap_ach/totPotential * 100;


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


end

