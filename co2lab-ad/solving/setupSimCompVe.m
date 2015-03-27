function s = setupSimCompVe(Gt, rock2D, varargin)
% Wrapper for setupSimComp to produce a vertically - integrated version.
   
   opt.useNewStandard = false; 
   opt.trans = []; 
   opt = merge_options(opt, varargin{:}); 
   
   % Compute vertially-integrated transmissibilities if not provided
   if isempty(opt.trans)
      rock_tmp      = rock2D; 
      rock_tmp.perm = rock2D.perm .* Gt.cells.H; 
      T             = computeTrans(Gt, rock_tmp); 
      cf            = Gt.cells.faces(:, 1); 
      nf            = Gt.faces.num; 
      T             = 1 ./ accumarray(cf, 1 ./ T, [nf, 1]); 
   else
      T = opt.trans; 
   end
   
   T_all = T; 
   
   % Computing vertically-integrated pore - volume
   pv = poreVolume(Gt, rock2D); 
   pv = pv .* Gt.cells.H; 
   
   % Calling setupSimComp, substituting the pore-volume and transmissiblity
   % values it computes with our vertically-integrated quantities.
   if opt.useNewStandard
      s = setupOperatorsTPFA(Gt, rock2D, 'porv', pv, 'trans', T); 
   else
      warning(['Using old standard for the operator structure, which uses ' ...
             'different signs for gradient and divergence term.  Consider ' ...
             'migrating to the new standard, as implemented in ' ...
             '''setupOperatorsTPFA''']); 
       s = setupSimComp(Gt, rock2D, 'porv', pv, 'trans', T); 
   end

   s.T_all = T_all; 
end