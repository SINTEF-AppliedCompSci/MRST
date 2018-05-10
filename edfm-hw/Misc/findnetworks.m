function networks=findnetworks(fracplanes,varargin)
% Given fracplanes (provided preprocessed before), find networks of
% connected fractures. 'clusters' is a cell array containing sets of
% fracplane indices.

opt=struct('histogram',false);
opt=merge_options(opt,varargin{:});


assert(isfield(fracplanes,'intersects'),'EDFM preprocessing needs to be completed first');
networks=cell(0);

remaining=1:length(fracplanes); % initialize with all fracplane indices

while ~isempty(remaining)
    newnetwork=remaining(1);
    newmembers=remaining(1);
    while ~isempty(newmembers)
        % obtain fractures that are connected to new members in the cluster
        connectedfracs=vertcat(fracplanes(newmembers).intersects)'; 
        connectedfracs=unique(connectedfracs);
        
        % filter out connected fractures that are already in the cluster
        connectedfracs=setdiff(connectedfracs,newnetwork);
        
        % assign connected fractures as new members
        newmembers=connectedfracs;
        
        % add connected fractures to new cluster
        newnetwork=[newnetwork,newmembers];    
    end
    
    % save new cluster
    networks{end+1}=newnetwork;
    
    % remove cluster from set of remaining fractures
    remaining=setdiff(remaining,newnetwork);
end

%% Plotting Histogram
if opt.histogram
    clustersize=cellfun(@length,networks);
    clustersize=horzcat(clustersize);
    
    [binneddata,~]=histcounts(clustersize,0.5:length(fracplanes)+0.5);
%     binneddata=log10(binneddata);
    xaxis=(1:length(fracplanes))*100/length(fracplanes);
%     xaxis=log10(xaxis);
    figure;
%     xaxis=xaxis(binneddata>0);
%     binneddata=binneddata(binneddata>0);
    bar(xaxis,binneddata);
    axis tight;
    xlabel('Network % size (proportion of fractures)');
    ylabel('Number of networks');
    title('Histogram of Network Size');
    
end


end