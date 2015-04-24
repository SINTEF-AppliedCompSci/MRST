function fluid = addSampledFluidProperties(fluid, shortname, varargin)

   opt.pspan  = [0.1, 400] * mega * Pascal; % CO2 default pressure range
   opt.tspan  = [  4, 250] + 274;           % CO2 default temperature range
   opt.pnum   = 800; % number of pressure samples
   opt.tnum   = 800; % number of temperature samples
   opt.props  = [true false false]; % which props to include [rho, mu, h]
   opt.fixedT = [];
   opt.assert_range = false;
   opt.nan_outside_range = false;

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
      fluid.(['rho', shortname]) = load_property(opt, 'D', fluidname, opt.fixedT, ...
                                                 opt.assert_range, opt.nan_outside_range); 
   end
   if opt.props(2)
      fluid.(['mu' , shortname]) = load_property(opt, 'V', fluidname, opt.fixedT, ...
                                                 opt.assert_range, opt.nan_outside_range); 
   end
   if opt.props(3)
      fluid.(['h'  , shortname]) = load_property(opt, 'H', fluidname, opt.fixedT, ...
                                                 opt.assert_range, opt.nan_outside_range); 
      
      if opt.props(1) % we have both enthalpy and density - we can also
                      % include internal energy
         fluid.(['u', shortname]) = ...
             @(P, T) fluid.(['h', shortname])(P, T) - P./fluid.(['rho',shortname])(P, T);
      end
   end
end

% ----------------------------------------------------------------------------

function pfun = load_property(opt, pname, fluidname, fixedT, assert_range, nan_outside)

   tabledir = [fileparts(mfilename('fullpath')) '/sampled_tables/'];
   fname = [tabledir, propFilename(opt.pspan, opt.tspan, opt.pnum, opt.tnum, fluidname, pname)];
   
   if (exist(fname) ~= 2)
      % data table not yet generated.  The following command will generate
      % and save them.   @@ User need 'coolprops' for this to work!
      try
         generatePropsTable(tabledir, fluidname, pname, opt.pspan, opt.tspan, opt.pnum, ...
                            opt.tnum);
      catch ME
         if strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            error(['Failed to generate property tables, as CoolProps could not ' ...
                   'be found.  Make sure you have CoolProps installed with ' ...
                   'the Matlab wrapper, and that the directory of the ' ...
                   'corresponding ''PropsSI'' is in your Matlab search ' ...
                   'path.']);
         else
            error(['Generating new property table failed due to an unknown ' ...
                   'error.']);
         end
      end
   end

   % We here load the generated table, and construct an object with evaluator
   % functions 
   obj = SampledProp2D(pname, fname, ...
                       'assert_in_range', assert_range, ...
                       'nan_outside_range', nan_outside);
   
   % We return the main evaluator function (which also works in an ADI-setting)
   pfun = obj.([pname]);
   
   if ~isempty(fixedT)
      % Temperature should be considered fixed -> property becomes function
      % of pressure only.
      pfun = @(p) pfun(p, fixedT);
   end
   
end

