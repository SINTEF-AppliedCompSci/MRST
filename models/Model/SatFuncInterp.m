function [krw, kro, pc] = SatFuncInterp(T)
   check_table(T);
   krw = @(s) interp_relperm(s, T(:, [1, 2]));
   kro = @(s) interp_relperm(s, T(:, [1, 3]));
   if size(T, 2) == 4
      pc = @(s) interp_pcap(s, T(:, [1, 4]));
   else
      pc = @zero_capillary;
   end
end
%--------------------------------------------------------------------------
function varargout = interp_relperm(s, T)
%    sw = max(s(:,1), T(1,1));
   [varargout{1 : nargout}] = interpolate(T(:, 1),T(:, 2),value(s));
%    if nargout > 1,
%       varargout{2} = [ varargout{2}(:, 1)    , ...
%                       zeros([size(sw, 1), 2]), ...
%                       -varargout{2}(:, 2) ];
%    end
end
%--------------------------------------------------------------------------
function varargout = interp_pcap(s, T)
%     sw = max(s(:,1), T(1,1));
    [varargout{1 : nargout}] = interpolate(T(:, 1),T(:, 2),value(s));                                      
end
%--------------------------------------------------------------------------
function varargout = zero_capillary(s, varargin)
   varargout(1 : nargout) = { zeros([size(s, 1), 1]) };
end

%--------------------------------------------------------------------------
function varargout = interpolate(s, y, sw)
    varargout{1} = interp1(s, y, sw, 'linear');
%     [b, b] = histc(sw, [-inf; s(2 : end-1); inf]);    
%     if nargout > 1,
%         % Caller requested derivatives too.
%         DY = bsxfun(@rdivide, diff(y), diff(s));
%         varargout{2} = DY(b, :);
%     end
end
%--------------------------------------------------------------------------
function check_table(T)
   % Basic validation of input data.
   assert (~any(any(T(:, 1:3) < 0)), ...
           'Saturation and rel-perm values must be non-negative.');
   assert (~ (T( 1 ,2) > 0), ...
           'Water must be immobile at connate water saturation.');
   assert (~ (T(end,3) > 0), ...
           'Oil must be immobile at maximum water saturation.');
   assert (all (diff(T(:,1)) > 0), ...
           'Water saturation must be monotonically increasing.');
   assert (~any(diff(T(:,2)) < 0), ...
           'Water rel-perm must be level or increasing down column.');
   assert (~any(diff(T(:,3)) > 0), ...
           'Oil rel-perm must be level or decreasing down column.');
end
%--------------------------------------------------------------------------

