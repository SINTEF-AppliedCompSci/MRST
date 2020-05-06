function [T, T_noflow] = computeMultiPointTrans(g, rock, varargin)

   opt = struct('useTensorAssembly', false, ...
                'neumann', true);
   
   [opt, extra] = merge_options(opt, varargin{:});
   
   if ~opt.useTensorAssembly
       [T, T_noflow] = computeMultiPointTransLegacy(g, rock, extra{:});
   else
       if opt.neumann
           mpfastruct = computeNeumannMultiPointTrans(G, rock, extra{:})
           T = mpfastruct;
       else
           mpfatstruct = computeMultiPointTransTensorAssembly(G, rock, extra{:});
           T = mpfastruct;
       end
   end       

end
