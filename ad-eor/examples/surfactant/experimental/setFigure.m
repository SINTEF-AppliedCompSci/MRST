function setFigure(fignum)
    try
        set(0, 'CurrentFigure', fignum);
    catch
        figure(fignum);
    end
end