function simpleUIbar(parent, data, start, height, txt)
    h = height/2;
    set(0, 'CurrentFigure', parent);
    axes('position', [0.1, start + h, .8, h])
    plotProgress(data)
    
    uicontrol(parent, 'units', 'Normalized', ...
        'Style', 'text', 'string', txt,...
        'position', [0.1, start, .8, h*0.9])
end

function plotProgress(data)
    rectangle('Position', [0, 0.01, data, 0.99], 'FaceColor', 'r')
    hold on
    rectangle('Position', [0, 0.01, 1, 0.99], 'FaceColor', 'none')
    ylim([0, 1]);
    xlim([0, 1]);
    axis off
    txt = sprintf('%2.1f%%', 100*data);
    text(0.5, 0.5, txt)
end
