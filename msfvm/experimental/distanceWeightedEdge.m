function newperm = distanceWeightedEdge(G, edge, perm, varargin)
% Alters permeability for some subset of nodes for some metric in some
% space
%
% G: Grid
% edge: The edge cells designated for change
% perm: permeability (or some other cell wise value)
%
% Optional:
% distance: Either a value or a vector which signifies the distance in
% which nodes will be sampled
% norm: Some norm capable of operating on matrices
% space: logical for logical, else actual coordinates
   opt = struct(...
       'distance', [3 3 1],...
       'norm', @(val) sqrt(sum(val.^2, 2)), ...
       'verbose', mrstVerbose,...
       'space', 'logical'...
   );
   opt = merge_options(opt, varargin{:});


   opt.distance = repmat(opt.distance, 1, G.griddim - numel(opt.distance) + 1);


   distall = repmat(opt.distance, G.cells.num, 1);

   Nc = G.cells.num;
   if strcmpi(opt.space, 'logical')
        [i j k] = ind2sub(G.cartDims, 1:Nc);
        coords = [i' j' k'];
        coords = coords(:,1:G.griddim);
   else
        coords = G.cells.centroids;
   end

   newperm = perm;

   for i = 1:numel(edge)
       dispif(opt.verbose, sprintf('%d of %d\n', i, numel(edge)));
       c = edge(i);
       dist = abs(coords - repmat(coords(c,:), Nc, 1));
%         dist = (coords - repmat(coords(c,:), Nc, 1));
       epsilon = all(dist<=distall,2);
       dist_eps = dist(epsilon,opt.distance~=0)./repmat(opt.distance(opt.distance~=0), sum(epsilon),1);
       norm_eps = opt.norm(dist_eps);
       weights = 1 - norm_eps./max(norm_eps);
%        newperm(c,:) = sum(perm(epsilon,:).*repmat(weights, 1, G.griddim))/sum(weights);
        vals = repmat(weights, 1, G.griddim)./perm(epsilon,:);
       newperm(c,:) = sum(weights)./sum(vals);

   end
end
