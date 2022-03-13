function u = ifcond(u,v,cond) 
% this function should be expanded
%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

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
if(isa(u,'ADI') || isa(v,'ADI'))
    if(isa(u,'ADI') && isa(v,'ADI'))
        value=ifcond(u.val,v.val,cond);
        %inx=cond+1;
        u  = ADI(value, plusJac(lMultDiag(double(cond), u.jac),...
            lMultDiag(double(~cond), v.jac)));
    elseif ~isa(u,'ADI')  %u is a vector
        %[value, inx] = max([u v.val], [], 2);
        value=ifcond(u,v.val,cond);
        %inx=cond+1;
        u  = ADI(value, lMultDiag(double(~cond), v.jac));
    elseif ~isa(v,'ADI') %v is a vector
        u = ifcond(v,u,~cond);
    else        
        error('Not yet implemented ...');
    end
else    
    assert(all(size(v)==size(u)));
    if(any(~cond))
        u(~cond)=v(~cond);
    end
    %error('Not yet implemented ...');
end
end
function J = lMultDiag(d, J1)
n = numel(d);
D = sparse((1:n)', (1:n)', d, n, n);
J = cell(1, numel(J1));
for k = 1:numel(J)
    J{k} = D*J1{k};
end
end
function J = plusJac(J1, J2)
nv1 = size(J1{1},1);
nv2 = size(J2{1},1);
if  nv1 == nv2
    J = cellfun(@plus, J1, J2, 'UniformOutput', false);
else     % only other legal option is that nv1 = 1 or nv2 =1
    if nv1 == 1
        J = cell(1, numel(J1));
        for k = 1:numel(J)
            J{k} = repmat(J1{k}, [nv2, 1]) + J2{k};
        end
    else % nv2 = 1
        J = plusJac(J2, J1);
    end
end
end