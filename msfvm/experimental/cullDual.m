function dg = cullDual(g, dg, A)
    inner = g.faces.neighbors(:,1) ~= 0 & g.faces.neighbors(:,2) ~= 0;
    tmp = double(g.faces.neighbors(inner,:));
    A = sparse(vertcat(tmp(:,1), tmp(:,2), (1:g.cells.num)'),...
                vertcat(tmp(:,2), tmp(:,1), (1:g.cells.num)'),1);

   cn = g.cells.num;

   inner        = false(cn,1);
   inner(dg.ii) = true;

   A = double(A | (A*A == 2));


   no_inner = numel(dg.ii);
   while 1
%        new = A*inner;
%        inner

       A_inner = A(inner, inner);
       comp_inner = components(A_inner);
       comp_unique = unique(comp_inner);
       nci = numel(comp_unique);
       fprintf('Inner comp: %d \n', nci);
       tmp = false(cn, nci);
       inner_ind = find(inner);
%        clf;
       for i = 1:nci
           ind = false(cn, 1);
           ind(inner_ind(comp_inner == comp_unique(i))) = true;
%            plotGrid(g, ind, 'facec', 'blue');
           extend =  A*(ind) > 0;% & ~ind;
           tmp(:, i) = extend;
%            comp_inner(extend) = comp_unique(i);
       end
       tmp2 = A*tmp & ~tmp;
       t2 = sum(tmp2, 2);
       t = sum(tmp,2);

       inner = inner | t == 1 & ~(t & t2);

%        inner = inner | t == 1 & ~(t & t2);

%        clf;
%        plotGrid(g, inner, 'facea', .3);
%         plotCellData(g, components(A(inner, inner)), inner);
%        plotGrid(g, dg.lineedge, 'facec', 'red', 'facea', .3)
%        pause
        new_no_inner = sum(inner);
        if new_no_inner == no_inner
            break
        else
            no_inner = new_no_inner;
        end
   end
dg.ii = find(inner);
dg.ee = setdiff(1:g.cells.num, [dg.ii(:); dg.nn(:)])';

dg.lineedge = setdiff(dg.lineedge, dg.ii);
dg.ii  = setdiff(dg.ii, dg.nn);


end
% function bad = offLimits(A, inner, components, cn, Cn)
%     c = unique(components);
%
%     cnts = false(cn, Cn);
%     for i = 1:numel(c)
%         tmp = false(cn, 1);
%
%     end
%
% end
