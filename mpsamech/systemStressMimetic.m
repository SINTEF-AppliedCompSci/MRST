function s = systemStressMimetic(g, CC, varargin)

      
       %s = computeMimeticIP(g,rock,varargin{:});
       dims=g.griddim;
       no2dof =@(no) reshape(bsxfun(@plus,g.griddim*(no-1),[1:g.griddim])',[],1);
       
       s = computeStressMimeticIP(g,CC,'type','hybrid');
       stmp= computeStressMimeticIP(g,CC,'type','mixed');
       s.B=stmp.B;
       cellNo  = rldecode(1:g.cells.num, diff(g.cells.facePos), 2) .';
       s.C     = sparse(no2dof(1:numel(cellNo)), no2dof(cellNo), 1);
       s.D     = sparse(no2dof(1:numel(cellNo)), no2dof(double(g.cells.faces(:,1))), 1, ...
           numel(cellNo)*dims, g.faces.num*dims);
       s.sizeB = repmat(size(g.cells.faces, 1), [1,2]);
       s.sizeC = size(s.C);
       s.sizeD = size(s.D);

      
   
end


