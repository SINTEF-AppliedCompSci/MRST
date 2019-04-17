function rock = makeUniformRock(G,k,ra)
%Creat valid rock permeability that is uniform;
%   G - Grid structure of MRST
%   k - principal value of permeability tensor
%   ra - rotational angle of principal permeability tensor

assert((numel(k)==2&&numel(ra)==1)||(numel(k)==3&&numel(ra)==3));
if(numel(k)==2)
    ra=ra./180.*pi;
    R=[cos(ra) -sin(ra);sin(ra) cos(ra)];
    K=R*[k(1) 0;0 k(2)]*R';
    rock.perm=repmat([K(1,1) K(1,2) K(2,2)],G.cells.num,1);
else
    K=rotx(ra(1))*roty(ra(2))*rotz(ra(3))*...
        diag([k(1) k(2) k(3)])*rotz(ra(3))'*roty(ra(2))'*rotx(ra(1))';
    rock.perm=repmat([K(1,1) K(1,2) K(1,3) K(2,2) K(2,3) K(3,3)],G.cells.num,1);
end
end

