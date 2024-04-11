function model = setrhoGfun(model, varargin)
    
    value = varargin{end};
    
    % The following functions in the fluid object must change:
    % - rhoGS
    % - rhoG
    % - bG

    % we need to keep track of the original functions in order to 
    % construct modified versions.  If this is the first time we apply
    % this function, make new fields in the fluid object for the 
    % original functions/values.
    if ~isfield(model.fluid, 'rhoG_orig')
        model.fluid.rhoG_orig = model.fluid.rhoG;
        model.fluid.bG_orig = model.fluid.bG;
        model.fluid.rhoGS_orig = model.fluid.rhoGS;
    end
    
    % Create modified functions/values
    model.fluid.rhoG = @(p, varargin) value .* model.fluid.rhoG_orig(p, varargin{:});
    model.fluid.bG = @(p, varargin) value .* model.fluid.bG_orig(p, varargin{:});
    model.fluid.rhoGS = model.fluid.rhoGS_orig * value;
    
end
