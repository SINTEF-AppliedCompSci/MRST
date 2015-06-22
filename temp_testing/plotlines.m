function varargout = plotlines(fl,conn)
if isempty(conn), conn = size(fl,1); end
for i = 1:numel(conn)
x = fl(conn(i),1:2:3);
y = fl(conn(i),2:2:4);
line(x,y,'Color','k','LineWidth',1);  hold on
% text(mean(x),mean(y),num2str(conn(i)));
end
varargout{i} = gcf;
return