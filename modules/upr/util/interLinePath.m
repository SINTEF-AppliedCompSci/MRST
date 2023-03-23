function [p] = interLinePath(line, fh, lineDist, sePtn, interpol, varargin)
    % Interpolate a line path. 
    % Arguments:
    %   line        Coordinates of the fault line. Must be ordered.
    %   fh          A function handle for the distance function for the
    %               interpolation fh = 1 will give equividistant
    %               interpolation
    %   lineDist    Scalar which sets the initial guess for the distance
    %               between interpolation points (Relative to fh = 1)
    % varargin:    
    %               Arguments passed to fh
%{
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) 2016 Runar Lie Berge. See COPYRIGHT.TXT for details.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%}
    % Parameters
    TOL = 1e-4; maxIt = 10000;
    if size(line,1)==1
        p = line;
        return
    elseif isempty(line)
        p = [];
        return
    end
    if interpol
        % Create initial points, equally distributed.
        p       = eqInterpret(line, lineDist, sePtn);
        flag    = false(size(p,1),1);
        ptsUsed = [1; size(line,1)];
    else
        [p, origIx, ptsUsed] = subdivideLineSegments(line, lineDist, sePtn);
        flag = false(size(p,1),1); flag(origIx) = true;
        ptsUsed([1 end]) = [1 size(line,1)];
    end
    
    % add auxillary points
    if sePtn(1)~=0 
        p = [line(1,:);p];  flag=[true; flag]; end
    if sePtn(2)~=0
        p = [p;line(end,:)]; flag=[flag; true]; end
 
    % plot(line(:,1),line(:,2),':k','LineWidth',1), hold all,              % makePlot
    % ncol = size(p,1);                                                    % makePlot
    % col = tatarizeMap(ncol);                                             % makePlot

    count=0;
    while count<maxIt
      count = count+1;
      % Calculate distances, and wanted distances
      d    = distAlLine(line, p);
      pmid = (p(1:end-1,:) + p(2:end,:))/2;
      dw   = fh(pmid,varargin{:});
      if sePtn(1)~=0, dw(1)   = dw(1).*sePtn(1);   end
      if sePtn(2)~=0, dw(end) = dw(end).*sePtn(2); end

      % Determine whether we need to insert or remove points in each of the
      % line segments
      dd = d-dw;
      ind = [0; find(flag(2:end-1)); numel(d)];
      insert = zeros(numel(ind)-1,1);
      for i=1:numel(insert)
          ix   = ind(i)+1:ind(i+1);
          dsum = sum(dd(ix)); 
          if dsum>min(dw(ix))
              [~, id]   = max(dd(ix));
              insert(i) = id + ind(i);
          elseif dsum<-max(dw(ix))
              [~, id]   = min(dd(ix));
              id        = id + ind(i) + (id==1);
              insert(i) = -id * (~flag(id));
          end
      end
      % Insert or remove the points determined in the previous for-loop
      if any(insert~=0)
          offset = 0;
          for i=1:numel(insert)
              if insert(i)>0
                  id     = insert(i)+offset;
                  p      = [p(1:id,:); pmid(insert(i),:); p(id+1:end,:)];
                  flag   = [flag(1:id); false; flag(id+1:end)];
                  offset = offset+1;
                  % ncol = ncol+1;                                         % makePlot
                  % tmp = tatarizeMap(ncol);                               % makePlot
                  % col = [col(1:id,:); tmp(ncol,:); col(id+1:end,:)];     % makePlot
              elseif insert(i)<0
                  id     = -insert(i)+offset;
                  p      = p([1:id-1,id+1:end],:);
                  flag   = flag([1:id-1,id+1:end]);
                  offset = offset-1;
                  % col = col([1:id-1,id+1:end],:);                        % makePlot
              end
          end
          continue
      end
      
      % If we only have external nodes, we can do nothing.
      if size(p,1)<=2, return, end
      
      % for i=1:size(p,1)                                                  % makePlot
      %     plot(p(i,1),p(i,2),'.','MarkerSize',18,'color',col(i,:));      % makePlot
      % end                                                                % makePlot
      
      % Move points based on desired length
      Fb = dw - d;                       % Bar forces
      Fn = Fb(1:end-1) - Fb(2:end);      % Force on internal nodes
      moveNode = Fn*0.2;                 % Movement of each internal node.
      d = d + [moveNode(1); moveNode(2:end) - moveNode(1:end-1); -moveNode(end)];
      p = interpLine(line, d, find(flag), ptsUsed); % Update node positions

      % Terminate if Nodes have moved (relative) less  than TOL
      if all(abs(moveNode(~flag(2:end-1)))<TOL*lineDist), break; end
    end
    % for i=1:size(p,1)                                                    % makePlot
    %     plot(p(i,1),p(i,2),'.','MarkerSize',30,'color',col(i,:));        % makePlot
    % end                                                                  % makePlot
    
    if sePtn(1)~=0, p = p(2:end,:);end
    if sePtn(2)~=0, p = p(1:end-1,:);end
    if count == maxIt
        warning('Fault interpolation did not converge.')
    end

end


function [d] = distAlLine(line, p)
    % Calculates the distace between consecutive interpolation points along
    % line
    % Arguments:
    %   line    line that is interpolated
    %   p       Interpolation points
    % Returns:
    %   d       distance between consecutive points of p, along line

    TOL = 50*eps;
    
    N = size(p,1);
    d = zeros(N-1,1);
    jointDist = 0;
    for i = 1:size(line,1)-1
        lineStart = repmat(line(i,:), size(p,1),1);
        lineEnd = repmat(line(i+1,:), size(p,1),1);
        distA = eucDist(lineStart, p) + eucDist(p,lineEnd);
        distB = eucDist(lineStart,lineEnd);
        indx  = find(abs(distA - distB) < TOL*distB);
        if numel(indx)==0 
            jointDist = jointDist + eucDist(line(i,:), line(i+1,:));
            continue
        elseif numel(indx)>=2
            d(indx(1:end-1)) = sqrt(sum((p(indx(1:end-1),:) ... 
                             - p(indx(2:end),:)).^2,2));
            
        end
        if indx(1)>1 && eucDist(line(i,:),p(indx(1),:))>TOL
            d(indx(1)-1) = jointDist + eucDist(line(i,:), p(indx(1),:));
        end
        jointDist = eucDist(p(indx(end),:), line(i+1,:));
    end
end


function [d] = eucDist(a, b)
    d = sqrt(sum((a - b).^2,2));
end


function [newPoints] = interpLine(path, dt, fix, ptsUsed)
    distS = sqrt(sum(diff(path,[],1).^2,2));
    t = [0; cumsum(distS)];
    
    newPtsEval = [0; cumsum(dt)];
    if numel(fix)>2
        newPtsEval(fix) = t(ptsUsed);  % Points that cannot move
    else
        newPtsEval(end) = t(end);
    end
    newPoints = interp1(t,path,newPtsEval);
end