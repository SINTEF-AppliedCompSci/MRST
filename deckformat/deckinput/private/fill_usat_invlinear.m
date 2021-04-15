function tbl = fill_usat_invlinear(tbl)
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

   nrec = diff(tbl.pos);
   i    = nrec == 1;

   if any(i),
      src  = find(~i);
      dst  = find( i);

      [src_ix, ncpy] = rlencode(sum(bsxfun(@gt, ...
                                           dst, ...
                                           reshape(src, 1, [])), 2));
      src_ix         = src_ix + 1;

      nrow = nrec(src(src_ix)) - 1;

      [IA, JA] = blockDiagIndex(nrow);
      p0 = IA == JA + 0;  % Main diag
      pm = IA == JA + 1;  % First sub-diag

      ii = [IA(p0); IA(pm)];
      jj = [JA(p0); JA(pm)];
      vv = repmat(-1, [sum(pm), 1]);

      IB = rldecode(cumsum([1; nrow(1 : end - 1)]), ncpy);
      JB = 1 : sum(ncpy);

      [IX, JX] = blockDiagIndex(nrow, ncpy);
      iX = sub2ind([sum(nrow), sum(ncpy)], IX, JX);

      i1 = mcolon( tbl.pos(src(src_ix) + 0), ...
                  (tbl.pos(src(src_ix) + 1) - 1) - 1);
      i2 = i1 + 1;

      t1 = tbl.data(i1, :);  % 1:end-1
      t2 = tbl.data(i2, :);  % 2:end

      delta = sparse(IA(IA >= JA), JA(IA >= JA), 1) * (t2(:,1) - t1(:,1));
      alpha = cumprod(t1(:, 2:end), 2) ./ ...
              cumprod(t2(:, 2:end), 2);                        % {1\over 2}
      rhs   = cumprod(tbl.data(tbl.pos(dst), 2:end), 2);

      derived = zeros([dot(ncpy, nrow), size(rhs, 2)]);
      for var = 1 : size(rhs, 2),
         A = sparse(ii, jj, [alpha(:,var); vv]);
         B = sparse(IB, JB,  rhs  (:,var),       sum(nrow), JB(end));

         X = A \ B;

         derived(:, var) = X(iX);
      end

      x       = tbl.data(tbl.pos(dst(JX)), 1) + delta(IX);
      derived = [x, derived(:,1), derived(:,2) ./ derived(:,1)];

      nrec2      = nrec;
      nrec2(dst) = rldecode(nrec(src(src_ix)), ncpy);
      p2         = cumsum([1; nrec2]);

      data2 = zeros([p2(end) - 1, size(tbl.data, 2)]);
      i1    = mcolon(p2(1 : end - 1)    , p2(1 : end - 1) + nrec - 1);
      i2    = mcolon(p2(dst) + nrec(dst), p2(dst + 1)            - 1);

      data2(i1,:) = tbl.data;
      data2(i2,:) = derived;

      tbl.pos  = p2;
      tbl.data = data2;
   end
end
