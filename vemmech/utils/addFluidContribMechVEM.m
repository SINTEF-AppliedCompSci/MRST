function fbc = addFluidContribMechVEM(G, bc, rock, isdirdofs)
   if(isfield(bc, 'type'))
      assert(all(strcmp(bc.type, 'pressure')));
   end
   if(~isempty(bc))
      %normal = model.G.faces.normals(bc.face, :);
      %sgn = 2*(model.G.faces.neigbours(bc.face, :) == 0)-1;% normal should be in to domain
      %bcf = bc_force;
      % find weights
      faces    = bc.face;
      lnn      = G.faces.nodePos(faces+1)-G.faces.nodePos(faces);
      inodes   = mcolon(G.faces.nodePos(faces), G.faces.nodePos(faces+1)-1);
      nodes    = G.faces.nodes(inodes);
      lfacenum = rldecode([1:numel(faces)]', lnn);%#ok
      facenum  = rldecode(faces, lnn);
      
      %assert boundary face
      assert(all(sum(G.faces.neighbors(faces, :) == 0, 2) == 1));
      lcells = sum(G.faces.neighbors(faces, :), 2);
      N = bsxfun(@times, G.faces.normals(faces, :), (2*(G.faces.neighbors(faces, 1) == lcells)-1));

      if(G.griddim == 2)
         qf = N/2;
         qf = qf(lfacenum, :);
      else
         % here one should weith the normal with teh subface area
         % this could have been taken from precomputed weights??
         [qc qfs] = calculateQC_vec(G);%#ok
         N        = bsxfun(@rdivide, N, G.faces.areas(faces));
         qf       = bsxfun(@times, N(lfacenum, :), qfs(inodes));
         %error()
      end
      
      % find nodes corresponding to forces
      dofs    = mcolon(G.griddim*(nodes-1)+1, G.griddim*(nodes-1)+G.griddim);
      assert(all(any(G.faces.neighbors(facenum, :) == 0, 2)));
      lfalpha = rock.alpha(sum(G.faces.neighbors(facenum, :), 2));%scale the face with alpha
      force   = bsxfun(@times, qf, lfalpha.*rldecode(bc.value, lnn));
      % transform to dofs numbering
      %assert(all(isdirdofs(dofs) == 0));
      force = force';
      ndof  = G.nodes.num*G.griddim;
      fbc   = accumarray(dofs', force(:)', [ndof, 1]);
      % map to current degees of freedom
      fbc   = fbc(~isdirdofs);
   else
      ndof  = G.nodes.num*G.griddim;
      fbc   = zeros(ndof, 1);
      fbc   = fbc(~isdirdofs);
   end
end
