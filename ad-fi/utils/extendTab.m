function T = extendTab(T, pm)
   if nargin < 2, pm = 1; end

   if ~ iscell(T),
      T = extend(T, pm);
   else
      T = cellfun(@(t) extend(t, pm), T, 'UniformOutput', false);
   end
end

%--------------------------------------------------------------------------

function T = extend(T, pm)
   t1 = T( 1 ,:); t1(1) = t1(1) - pm;
   te = T(end,:); te(1) = te(1) + pm;

   T = [t1; T; te];
end
