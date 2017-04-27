function plotFaces2D(G, faces, varargin)

opt=struct('col','r','uu',[],'shift',[],'linewidth',4);
opt=merge_options(opt,varargin{:});
if(~isempty(opt.uu))
    G.nodes.coords=G.nodes.coords+opt.uu;
end
opt=merge_options(opt,varargin{:});
if(G.griddim==2)
      nodes=G.faces.nodes(mcolon(G.faces.nodePos(faces),G.faces.nodePos(faces+1)-1));
      nodes=reshape(nodes,2,[])';
      
      for j=1:size(nodes,1)
          %line(G.nodes.coords(nodes(1,:),1),G.nodes.coords(nodes(j,:),2),varargin{:})
          pos1=G.nodes.coords(nodes(j,:),1);
          pos2=G.nodes.coords(nodes(j,:),2);
          if(~isempty(opt.shift))
          pos1=pos1+bsxfun(@times,G.faces.normals(faces(j),1),opt.shift(j));
          pos2=pos2+bsxfun(@times,G.faces.normals(faces(j),2),opt.shift(j));
          end
          line(pos1,pos2,'Color',opt.col,'LineWidth',opt.linewidth);
      end
      return
else
    error('This function do only work for 2D')
end
