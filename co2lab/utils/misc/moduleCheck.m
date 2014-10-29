function moduleCheck(varargin)
%% Load modules whose names are in the argument list, unless they are loaded already.

   for module = reshape(varargin, 1, []),
      m = module{1};

      try
         require(m);
      catch
         fprintf('Loading module %s\n', m);
         mrstModule('add', m);
      end
   end
end
