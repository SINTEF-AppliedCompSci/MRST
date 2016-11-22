function bc = addBCVEM(bc, f, t, g, varargin)

opt = struct('sat', []);
opt = merge_options(opt, varargin{:});
s   = opt.sat;

if isempty(bc),
   bc = struct('face'  , {{}}  , ...
               'type' , {{}}, ...
               'func', {{}}  , ...
               'sat'   , []);
end

nn = numel(bc.face)+1;
bc.face{nn} = f;
bc.type{nn} = t;
bc.func{nn} = g;

end