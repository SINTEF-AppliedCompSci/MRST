function [aq, unitNo] = getRestartAquiInfo(rstrt, step)
%Undocumented Utility Function

%{
Copyright 2009-2019 SINTEF Digital, Mathematics & Cybernetics.

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

% currently only considers analytical aquifers
intehead = rstrt.INTEHEAD{step};
unitNo = intehead( 3);
aq = [];

naq = intehead(41);    % number of aquifers (guess)
if naq > 0 && isfield(rstrt, 'XAAQ')
    nxaaq = intehead(45);
    xaaq  = rstrt.XAAQ{step};
    if ~(numel(xaaq) == naq*nxaaq)
        warning('Unexpected format of aquifer data in restart file, skipping ...')
        return;
    else
        aq = repmat(struct('cijk', [], 'p', [], 'qW', [], 'flux', [], 'num', []), [naq, 1]);
        for k = 1:naq
            aq(k).p  = xaaq(2 + (k-1)*nxaaq);
            aq(k).qW = xaaq(1 + (k-1)*nxaaq);
        end
        % get number of connections for each aquifer
        niaaq = intehead(43);
        nconn = rstrt.IAAQ{step}(1 + (0:(naq-1))*niaaq);
        
        nicaq = intehead(46);
        nacaq = intehead(48);
        
        % if naq>1, fields will have number suffix (_1, _2, ...)
        sfx = @(k)'';
        if naq > 1
            sfx = @(k)sprintf('_%d', k);
        end
        for k = 1:naq
            icaq = rstrt.(['ICAQ', sfx(k)]){step};
            dc   = (0:(nconn(k)-1))*nicaq;
            aq(k).cijk = icaq(bsxfun(@plus, dc(:), 1:3));
            
            acaq = rstrt.(['ACAQ', sfx(k)]){step};
            d = (0:(nconn(k)-1))*nacaq;
            aq(k).flux = acaq(1+d);
            
            acaqnum = rstrt.(['ACAQNUM', sfx(k)]){step};
            aq(k).num = acaqnum;
        end
    end
    
end
    

end
