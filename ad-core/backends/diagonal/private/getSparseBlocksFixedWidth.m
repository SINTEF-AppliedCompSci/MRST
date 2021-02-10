function [I, J, V, imax, jmax] = getSparseBlocksFixedWidth(D)
              % TODO: Check for row major version
            V = D.diagonal;
            if D.rowMajor
                V = V';
            end
            n = size(V, 1);
            m = D.dim(2);
            nmap = size(D.map, 2);
            imax = n;
            jmax = prod(D.dim);
            if isempty(D.subset)
                nval = n;
            else
                nval = numel(D.subset)/nmap;
            end
            if size(D.map, 1) == 1
                D.map = repmat(D.map, nval, 1);
            end
            I = repmat((1:n)', 1, nmap*m);
            jmap = D.map;
            if ~isempty(D.subset)
                jmap = jmap(D.subset, :);
            end
            if ~isempty(D.parentSubset)
                jmap = D.parentSubset(jmap);
            end
            
            J = jmap(:, rldecode((1:nmap)', m));
            tmp = (0:m-1)*D.dim(1);
            J = bsxfun(@plus, J, repmat(tmp, 1, nmap));
            V(all(jmap == 0, 2), :) = 0;
  end