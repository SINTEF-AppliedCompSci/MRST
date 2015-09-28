function [ hfig_new ] = setColorbarHandle( hfig, varargin )
% Function to ensure pre-R2014 and post-R2014 compatibility of handles to
% figures, colorbars, and axes

    

opt.LabelName   = 'defaultLabelName';
opt.fontSize    = 14;

opt = merge_options(opt, varargin{:});


figure(hfig);
hcb_curr = colorbar;


    % ---------------------------------------------------------------------
    % Colorbar handles:
    
    if isa(hcb_curr, 'matlab.graphics.illustration.ColorBar')
        % This is 2014b or later
        
        hcb_curr.Label.String = opt.LabelName;
        
        
    elseif isfloat(hcb_curr)
        % This is 2014a or earlier
        
        %set(hcb_curr, 'XTickLabel', opt.LabelName);
        ylabel(hcb_curr, opt.LabelName, 'fontSize', opt.fontSize);
        
    end

    
    set(hcb_curr, 'fontSize', opt.fontSize);
    

hfig_new = gcf;



end

