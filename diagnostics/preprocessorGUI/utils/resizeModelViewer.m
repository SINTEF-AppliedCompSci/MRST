function resizeModelViewer(src, event)
fig = src;
c = fig.Children;
histix = find(strcmp('histogram', get(c,'Tag')));
cbarix = find(strcmp('Colorbar', get(c,'Tag')));

for k = 1:numel(histix)
    c(histix(k)).Position(1) = sum(c(cbarix(k)).Position([1 3])) + 0.002;
end
end

