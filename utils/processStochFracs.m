function [set1,nonPlanarSets,fracArea] = processStochFracs(fracSet)
%REMOVENANS Summary of this function goes here
%   Detailed explanation goes here
    set1 = {};
    ct=0;
    jj=0;
    kk=0;
    nonPlanarSets = {};
    fracArea = zeros(numel(fracSet),1);
    for i=1:numel(fracSet)
        %remove nans from nodal coordinates
        temp = fracSet{i}(1:end-1,:);
        fracArea(i) = convexpolygonarea(temp,1e-6);
        
        if sum(isnan(temp))
            continue
        end
        
        %remove fractures that are too small!
        if (fracArea(i)<1e-5)
            continue
        end
        
        %remove lines, fracs should be at least 2D
        [lenOfTemp, ~] = size(temp);
        if lenOfTemp < 3
            continue
        end
            
        %remove non-planar fractures
        [fem, disit] = iscoplanar(temp,'distTolerance', 1e-9);
        if (all(fem) ) 
            ct = ct+1;
            set1{ct} = fracSet{i};
        else
            if(all(disit))
                kk=kk+1;
                ct = ct+1;
                set1{ct} = fracSet{i};
            else
                jj = jj + 1;
                nonPlanarSets{jj} = fracSet{i};
            end
        end
            
    end
    set1 = set1';
    nonPlanarSets=nonPlanarSets';
    fprintf('%d  out of %d were coplanar based on disit but were not coplanar based on Delaunay triangulation\n',kk,numel(fracSet));
end

