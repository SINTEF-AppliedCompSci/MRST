function fluidimp = veFluid2Fluid(G,fluid,H)   
% Function to convert ve-fluid (made with initVEFluid to fluid that can be
% used with s-formulation, at the moment hack and hardcoded geometry

%{
#COPYRIGHT#
%}

      fluidimp.pc=@(rsol) pc(rsol, fluid, H);
      fluidimp.relperm =@(s, varargin) relperm(s);
      
      fluidimp.saturation = @(rsol) rsol.s;  
      fluidimp.properties = @(varargin) properties(fluid,varargin{:});

end

function varargout = pc(rsol, fluid, H)   
   varargout{1}    =  norm(gravity)*(fluid.rho(1)-fluid.rho(2))*fluid.pc(rsol.s.*H);

   if nargout > 1,
      varargout{2} = norm(gravity)*(fluid.rho(1)-fluid.rho(2))*H.*ones(numel(rsol.s),1);
   end
end


function [mu, rho] = properties(fluid,varargin)
   mu = fluid.mu; 
   rho =fluid.rho;
end
function [kr, dkr] = relperm(s)
   kr = [s,1-s]; 
   dkr = repmat([1,-1],numel(s),1);
end
