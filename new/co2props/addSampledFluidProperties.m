function fluid = addSampledFluidProperties(fluid, shortname, varargin)

   opt.pspan  = [10 410] * barsa;
   opt.tspan  = [10 500] + 273.15;
   opt.pnum   = 200; % number of pressure samples
   opt.tnum   = 200; % number of temperature samples
   opt.ref    = [100 * barsa, 35 + 273.15]; % BW reference temperature/pressure
   opt.props  = [true false false]; % which props to include [rho, mu, h]
   opt.fixedT = [];

   opt = merge_options(opt, varargin{:});
   
   % Determine full fluid name 
   switch shortname
     case 'G'
       fluidname = 'CarbonDioxide';
     case 'W'
       fluidname = 'Water';
     otherwise
      error('Unsupported fluid name');
   end
   
   % Add density, viscosity and enthalpy properties
   if opt.props(1)
      fluid.(['rho', shortname]) = load_property(opt, 'D', fluidname, opt.fixedT); 
   end
   if opt.props(2)
      fluid.(['mu' , shortname]) = load_property(opt, 'V', fluidname, opt.fixedT); 
   end
   if opt.props(3)
      fluid.(['h'  , shortname]) = load_property(opt, 'H', fluidname, opt.fixedT);
      
      if opt.props(1) % we have both enthalpy and density - we can also
                      % include internal energy
         fluid.(['u', shortname]) = ...
             @(P, T) fluid.(['h', shortname])(P, T) - P./fluid.(['rho',shortname])(P, T);
      end
   end
end

% ----------------------------------------------------------------------------

function pfun = load_property(opt, pname, fluidname, fixedT)

   fname = propFilename(opt.pspan, opt.tspan, opt.pnum, opt.tnum, fluidname, pname);
   
   if (exist(fname) ~= 2)
      % data table not yet generated.  The following command will generate
      % and save them.   @@ User need 'coolprops' for this to work!
      generatePropsTable(fluidname, pname, opt.pspan, opt.tspan, opt.pnum, opt.tnum);
   end

   % We here load the generated table, and construct an object with evaluator
   % functions 
   obj = SampledProp2D(pname, fname);
   
   % We return the main evaluator function (which also works in an ADI-setting)
   pfun = obj.([pname]);
   
   if ~isempty(fixedT)
      % Temperature should be considered fixed -> property becomes function
      % of pressure only.
      pfun = @(p) pfun(p, fixedT);
   end
   
end

