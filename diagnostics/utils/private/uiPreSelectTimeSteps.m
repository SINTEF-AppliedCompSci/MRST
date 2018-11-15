function steps = uiPreSelectTimeSteps(info)
d = dialog('Resize', 'on');
s = PreDiagnosticsSelector('Parent', d, 'restartInfo', info);
d.CloseRequestFcn = @setSteps;
d.SizeChangedFcn  = @updateSize;
uiwait(d)

    function setSteps(src, event)
    steps = s.ix;
    delete(d)
    end

    function updateSize(src, event)
        s.Position = [0 0 d.Position(3:4)];
    end
end