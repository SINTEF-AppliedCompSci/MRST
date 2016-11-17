function [p] = interLinePath(line, fh, lineDist,sePtn, varargin)
    % Interpolate a line path. 
    % Arguments:
    %   line        Coordinates of the fault line. Must be ordered.
    %   fh          A function handle for the distance function 
    %               for the interpolation fh = 1 will give a equiv distant
    %               interpolation
    %   lineDist    Scalar which set the initial guess for the distance 
    %               between interpolation
    %               points (Relative to fh = 1)
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
    % Create initial points, equally distributed.
    p = eqInterpret(line, lineDist,sePtn);
    % add auxillary points
    if sePtn(1)~=0, p = [line(1,:);p];   end
    if sePtn(2)~=0, p = [p;line(end,:)]; end
    count=0;
    while count<maxIt
      count = count+1;
      % Calculate distances, and wanted distances
      d = distAlLine(line, p);
      pmid = (p(1:end-1,:) + p(2:end,:))/2;
      dw = fh(pmid,varargin{:});
      if sePtn(1)~=0, dw(1)   = dw(1).*sePtn(1);   end
      if sePtn(2)~=0, dw(end) = dw(end).*sePtn(2); end

      % Possible insert or remove points
      if sum(d - dw)>min(dw)
          [~, id] = max(d - dw);
          p = [p(1:id,:); pmid(id,:); p(id+1:end,:)];
          continue
      elseif sum(d - dw)<-max(dw)
          [~, id] = min(d - dw);
          if id == 1, id = 2; end
          p = p([1:id-1,id+1:end],:);
          continue
      end
      % If we only have external nodes, we can do nothing.
      if size(p,1)<=2, return, end
      % Move points based on desired length
      Fb = dw - d;                       % Bar forces
      Fn = Fb(1:end-1) - Fb(2:end);      % Force on internal nodes
      moveNode = Fn*0.2;                 % Movement of each internal node.
      d = d + [moveNode(1); moveNode(2:end) - moveNode(1:end-1); -moveNode(end)];
      p = interpLine(line,d);            % Update node positions

      % Terminate if Nodes have moved (relative) less  than TOL
      if all(abs(moveNode)<TOL*lineDist), break; end
    end
    
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


function [newPoints] = interpLine(path, dt)
    distS = sqrt(sum(diff(path,[],1).^2,2));
    t = [0; cumsum(distS)];
    
    newPtsEval = [0; cumsum(dt)];
    newPtsEval(end) = t(end); % Last point can not move

    newPoints = interp1(t,path,newPtsEval);
end