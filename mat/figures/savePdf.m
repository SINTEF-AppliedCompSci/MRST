function savePdf(h, fileName, varargin)

    opt = struct('cut' , 4);
    opt = merge_options(opt, varargin{:});
    cut = opt.cut;

    ps = get(h, 'Position');
    ratio = (ps(4)-ps(2)) / (ps(3)-ps(1));
    paperWidth = 10;
    paperHeight = paperWidth*ratio - cut;
    set(h, 'paperunits', 'centimeters');
    set(h, 'papersize', [paperWidth paperHeight]);
    set(h, 'PaperPosition', [0    0   paperWidth paperHeight]);

    print(h, '-dpdf', fileName);
    
end