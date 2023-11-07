function model = setrhoGmultfun(model, varargin)
    
    val = varargin{end};

    % The following functions in the fluid object must change:
    % - rhoG
    % - bG

    % At first call, we modify the fluid object to accomodate for the
    % changes
    if ~isfield(model.fluid, 'rhoGmult')
        % Keeping reference values of all fields that will be impacted by 
        % changes to 'rhoGmult'
        model.fluid.rhoGmult = 1;
        model.fluid.rhoGS_orig = model.fluid.rhoGS;
        model.fluid.rhoG_orig = model.fluid.rhoG;
        model.fluid.pcWG_orig = model.fluid.pcWG;
        model.fluid.bG_orig = model.fluid.bG;
    end
    
    % set multiplier and update density functions
    model.fluid.rhoGmult = val;
    model.fluid.rhoGS = model.fluid.rhoGS_orig .* val;
    model.fluid.rhoG = @(p, varargin) val .* model.fluid.rhoG_orig(p, varargin{:});
    
    % change in upscaled capillary pressure is best expressed through an
    % intermediary function
    cap_adjust = @(p) (model.fluid.rhoW(p) - val .* model.fluid.rhoG_orig(p)) ./ ...
                      (model.fluid.rhoW(p) - model.fluid.rhoG_orig(p));
    
    model.fluid.pcWG = @(sg, p, varargin) cap_adjust(p) .* ...
                                          model.fluid.pcWG_orig(sg, p, varargin{:});
    
    
end
