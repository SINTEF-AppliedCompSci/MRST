function varargout = plotFractureLines(G,fracture,varargin)
% This function plots the fracture lines using the fracture structure
% returned by processFracture2D

dc = zeros(numel(fracture.lines),3);
if ~isempty(varargin)
    if strcmp(varargin{1},'lines')
%       dc = distinguishable_colors(numel(fracture.lines),{'w'});
        dc = rand(numel(fracture.lines),3);
    elseif strcmp(varargin{1},'network')
%       dn = distinguishable_colors(numel(fracture.network),{'w'});
        dn = rand(numel(fracture.network),3);
        for i = 1:numel(fracture.network)
            lines = fracture.network(i).lines;
            dc(lines,:) = repmat(dn(i,:),numel(lines),1);
        end
    elseif isnumeric(varargin{1})
        dc = rand(numel(varargin{1}),3);
    end
end
h = plotGrid(G,'FaceColor','w','EdgeAlpha',0.07);
ax = gca;
hold on
if isnumeric(varargin{1})
    lines = varargin{1};
    for i = 1:numel(lines)
    x = [fracture.lines(lines(i)).endp(1);fracture.lines(lines(i)).endp(3)];
    y = [fracture.lines(lines(i)).endp(2);fracture.lines(lines(i)).endp(4)];
    line(x,y,'Color',dc(i,:),'Parent',ax,'LineWidth',1.5);
%     text(mean(x),mean(y),num2str(i));
    end
else
for i = 1:numel(fracture.lines)
    x = [fracture.lines(i).endp(1);fracture.lines(i).endp(3)];
    y = [fracture.lines(i).endp(2);fracture.lines(i).endp(4)];
    line(x,y,'Color',dc(i,:),'Parent',ax,'LineWidth',1.5);
%     text(mean(x),mean(y),num2str(i));
end
end
varargout{1} = h;
return 