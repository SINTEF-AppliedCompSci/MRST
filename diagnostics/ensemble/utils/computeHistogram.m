function h = computeHistogram(data,varargin)
% Histogram calculator
%
% SYNOPSIS:
%   histogram_structure = computeHistogram(data,varargin)
%
% DESCRIPTION:
%  
% REQUIRED PARAMETERS:
%   data - array that contains the data to obtains its histogram
%
% OPTIONAL PARAMETERS:
%   N_bins    -  number of bins for the histograms
%   Edges     -  is an array, specifies the edges of the bins.
%   LogEdges  -  is an array, specifies the edges of the bins for the 
%                logarithmic hiscogram.
%   BinLimits -  [BMIN,BMAX] bins only elements 
%                in X between BMIN and BMAX inclusive, 
%   Normalization pecifies the normalization scheme of the histogram values returned in N. NM can be:
%                    'count'   Each N value is the number of observations in 
%                              each bin. SUM(N) is generally equal to 
%                              NUMEL(X), but is less than if some of 
%                              the input data is not included in the bins. 
%                              This is the default.
%              'probability'   Each N value is the relative number of 
%                              observations (number of observations in bin / 
%                              total number of observations), and SUM(N) is 
%                              less than or equal to 1.
%             'countdensity'   Each N value is the number of observations in 
%                              each bin divided by the width of the bin. 
%                      'pdf'   Probability density function estimate. Each N 
%                              value is, (number of observations in bin) / 
%                              (total number of observations * width of bin).
%                 'cumcount'   Each N value is the cumulative number of 
%                              observations in each bin and all previous bins. 
%                              N(end) is less than or equal to NUMEL(X).
%                      'cdf'   Cumulative density function estimate. Each N 
%                              value is the cumulative relative number of 
%                              observations in each bin and all previous bins. 
%                              N(end) is less than or equal to 1.

opt = struct('N_bins', 500, ...
             'Edges',   [], ...
             'LogEdges',   [], ...
             'BinLimits',   [], ...
             'Normalization', 'cumcount');         
         
opt = merge_options(opt, varargin{:});

    minimun = min(data);

% Calculating max and min if they aren't input    
if isempty(opt.BinLimits)
    opt.BinLimits(1,:) = minimun;
    opt.BinLimits(2,:) = max(data);
end

    if isempty(opt.Edges)
        
        for i =  1: size(data,2)
             
             h.min(1,i) = opt.BinLimits(1,i);                          
             h.max(1,i) = opt.BinLimits(2,i);
             
             
               
            [h.hist.n_counts(:,i),h.hist.edges(:,i)] = histcounts(data(:,i),opt.N_bins,...
                                                                          'BinLimits',[h.min(1,i),h.max(1,i)],...
                                                                          'Normalization',opt.Normalization);

                                                                      
            % Check for minimun values to make it safe for logarithm
            if ((h.min(1,i)>0)&&(minimun(i)>0))
                [h.loghist.n_counts(:,i),edges] = histcounts(log10(data(:,i)),opt.N_bins,...
                                                                              'BinLimits',[log10(h.min(1,i)),log10(h.max(1,i))],...
                                                                              'Normalization',opt.Normalization);
                 h.loghist.edges(:,i) = 10.^edges;
                 
            elseif (h.min(1,i)==0)&&(minimun(i)>0)
                
                % remove zero elements for the logaritmic histogram
                 mini_aux = h.max(1,i)*10^(-5);
                 
                [h.loghist.n_counts(:,i),edges] = histcounts(log10(data(:,i)),opt.N_bins,...
                                                                              'BinLimits',[log10(mini_aux),log10(h.max(1,i))],...
                                                                              'Normalization',opt.Normalization);
                 h.loghist.edges(:,i) = 10.^edges;
            else
                 
                 i_log   = find(data(:,i)>0);
                 if ~isempty(i_log)
                    data_log_safe = data(i_log,i);
                    h.min(1,i) =  min(data_log_safe);
                    h.max(1,i) =  max(data_log_safe);
                    [h.loghist.n_counts(:,i),edges] = histcounts(log10(data_log_safe),opt.N_bins,...
                                                                                  'BinLimits',[log10(min(data_log_safe)),log10(max(data_log_safe))],...
                                                                                  'Normalization',opt.Normalization);
                     h.loghist.edges(:,i) = 10.^edges;
                 else
                     h.min(1,i) = 0;
                     h.max(1,i) = 1;
                     h.loghist.edges(:,i)    = h.hist.edges(:,i);
                     h.loghist.n_counts(:,i) = h.hist.n_counts(:,i);
                     warning('Logaritmic histogram is not available. There are zero values in data.');
                 end
            end
             
        end
    else 
        
        for i =  1: size(data,2)
            [h.hist.n_counts(:,i),h.hist.edges(:,i)] = histcounts(data(:,i),opt.Edges(:,i),...
                                                                    'Normalization',opt.Normalization);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
            if  ~isempty(opt.LogEdges)
                [h.loghist.n_counts(:,i),h.loghist.edges(:,i)] = histcounts(data,opt.LogEdges(:,i),...
                                                                    'Normalization',opt.Normalization);
            end
        end
    end

end


