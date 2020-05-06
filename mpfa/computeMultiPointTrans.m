function [T, T_noflow] = computeMultiPointTrans(G, rock, varargin)

   opt = struct('useTensorAssembly', false, ...
                'neumann', false);
   
   [opt, extra] = merge_options(opt, varargin{:});
   
   if ~opt.useTensorAssembly
       [T, T_noflow] = computeMultiPointTransLegacy(G, rock, extra{:});
   else
       if opt.neumann
           mpfastruct = computeNeumannMultiPointTrans(G, rock, extra{:})
           T = mpfastruct;
       else
           mpfastruct = computeMultiPointTransTensorAssembly(G, rock, extra{:});
           T = mpfastruct;
       end
   end       

end
