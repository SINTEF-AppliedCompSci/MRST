function fluid = makeADIFluid(type, use_residual, varargin)
   
   mu(3), rho(3), n(3), sr, sw, pvMultR, bW, bG, surface_tension

   %% Making the base fluid object that will be modified depending on type
   fluid = initSimpleADIFluid
      
      
   %% adding type-specific modifications
   switch type
     case 'simple' %@@ tested anywhere?
       fluid = makeSimpleFluid(fluid, use_residual);
     case 'integrated' %@@ tested anywhere? 
       fluid = makeIntegratedFluid(fluid, use_residual);
     case 'sharp interface'
       fluid = makeSharpInterfaceFluid(fluid, use_residual);
     case 'linear cap.'
       fluid = makeLinCapFluid(fluid, use_residual);     
     case 'S table'
       fluid = makeSTableFluid(fluid, use_residual);
     case 'P-scaled table'
       fluid = makePScaledFluid(fluid, use_residual);
     case 'P-K-scaled table'
       fluid = makePKScaledFluid(fluid, use_residual);
     otherwise
       error([type, ': no such fluid case.']);
   end
         
end

% ============================================================================

function fluid = makeSimpleFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

function fluid = makeIntegratedFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

function fluid = makeSharpInterfaceFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

function fluid = makeLinCapFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

function fluid = makeSTableFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

function fluid = makePScaledFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

function fluid = makePKScaledFluid(use_residual)
   
end

% ----------------------------------------------------------------------------

