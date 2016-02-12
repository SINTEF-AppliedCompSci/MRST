function Nc = computeCapillaryNumber(p, c, pBH, W, fluid, G, operators, varargin)

   opt = struct('velocCompMethod', 'linear');
   opt = merge_options(opt, varargin{:});

   s = operators;
   v = -s.T.*s.Grad(p);

   switch opt.velocCompMethod
     case 'linear'

       veloc = s.veloc;
       veloc_sq = 0;
       for i = 1 : numel(veloc)
          veloc_sq = veloc_sq + veloc{i}(v).^2;
       end
       [velocW, wc] = computeWellContrib(G, W, p, pBH);
       veloc_sq(wc) = veloc_sq(wc) + velocW.^2;

     case 'square'

       veloc_sq = s.sqVeloc(v);

       add_well_contrib = false;
       if add_well_contrib
          [velocW, wc] = computeWellContrib(G, W, p, pBH);
          veloc_sq(wc) = veloc_sq(wc) + velocW.^2;
       else

     otherwise
       error('option for velocCompMethod not recognized');

   end

   abs_veloc = (veloc_sq).^(1/2);

   % subplot(2, 1, 2)
   % plot(abs_veloc.val);
   % title('velocity')
   % drawnow;

   sigma = fluid.ift(c);
   Nc = abs_veloc./sigma;

end

function [dx, dy, dz] = cellDims(G, ix)
% cellDims -- Compute physical dimensions of all cells in single well
%
% SYNOPSIS:
%   [dx, dy, dz] = cellDims(G, ix)
%
% PARAMETERS:
%   G  - Grid data structure.
%   ix - Cells for which to compute the physical dimensions
%
% RETURNS:
%   dx, dy, dz -- [dx(k) dy(k)] is bounding box in xy-plane, while dz(k) =
%                 V(k)/dx(k)*dy(k)

    n = numel(ix);
    [dx, dy, dz] = deal(zeros([n, 1]));

    ixc = G.cells.facePos;
    ixf = G.faces.nodePos;

    for k = 1 : n,
       c = ix(k);                                     % Current cell
       f = G.cells.faces(ixc(c) : ixc(c + 1) - 1, 1); % Faces on cell
       e = mcolon(ixf(f), ixf(f + 1) - 1);            % Edges on cell

       nodes  = unique(G.faces.nodes(e, 1));          % Unique nodes...
       coords = G.nodes.coords(nodes,:);            % ... and coordinates

       % Compute bounding box
       m = min(coords);
       M = max(coords);

       % Size of bounding box
       dx(k) = M(1) - m(1);
       if size(G.nodes.coords, 2) > 1,
          dy(k) = M(2) - m(2);
       else
          dy(k) = 1;
       end

       if size(G.nodes.coords, 2) > 2,
          dz(k) = G.cells.volumes(ix(k))/(dx(k)*dy(k));
       else
          dz(k) = 0;
       end
    end
end

function [velocW, wc] = computeWellContrib(G, W, p, pBH)
% We use the bottom hole pressure to compute well velocity
   perf2well = getPerforationToWellMapping(W);
   Rw = sparse(perf2well, (1:numel(perf2well))', 1, numel(W), numel(perf2well));

   wc  = vertcat(W.cells);
   pW  = p(wc);
   pBH = Rw'*pBH;
   Tw  = vertcat(W.WI);
   Tw  = Rw'*Tw;
   welldir = { W.dir };
   i = cellfun(@(x)(numel(x)), welldir) == 1;
   welldir(i) = arrayfun(@(w) repmat(w.dir, [ numel(w.cells), 1 ]), ...
                         W(i), 'UniformOutput', false);
   welldir = vertcat(welldir{:});
   [dx, dy, dz] = cellDims(G, wc);
   thicknessWell = dz;
   thicknessWell(welldir == 'Y') = dy(welldir == 'Y');
   thicknessWell(welldir == 'X') = dx(welldir == 'X');

   if ~isfield(W, 'rR')
      error('The representative radius of the well is not initialized');
   end
   rR = vertcat(W.rR);

   velocW = Tw.*(pW - pBH)./(2 * pi * rR .* thicknessWell);

end