function [stress,strain]=calculateStressVEM(G,uu, op,varargin)
%Undocumented Utility Function

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

opt.do_patch=false;
opt=merge_options(opt,varargin{:});

if(G.griddim==2)
    lindim=3;
else
    lindim=6;
end
stress=reshape(op.D*op.WC'*op.assemb'*reshape(uu',[],1),lindim,[])';
if(opt.do_patch)
  stress=patchRecovery(G,stress);
end
if(G.griddim==3)
 stress=bsxfun(@rdivide,stress(:,[1:3,5,6,4]),[ones(1,3),2*ones(1,3)]);
 stress(:,[4,6]) = stress(:,[6,4]);
else
    assert(G.griddim==2)
  stress=bsxfun(@rdivide,stress,[ones(1,2),2]);  
end

if(nargout==2)
    error('Need to be checked')
    strain=reshape(op.WC'*op.assemb'*reshape(uu',[],1),lindim,[])';
end

end
