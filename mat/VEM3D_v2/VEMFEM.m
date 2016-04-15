
% % Q = orth(eye(NK)-PN);
% q1 = [ 1 -1 -1  1  1 -1 -1  1];
% q2 = [ 1 -1  1 -1 -1  1 -1  1];
% q3 = [ 1  1 -1 -1 -1 -1  1  1];
% q4 = [-1  1  1 -1  1 -1 -1  1];
% 
% hx = abs(max(X(:,1)) - min(X(:,1)))/2;
% hy = abs(max(X(:,2)) - min(X(:,2)))/2;
% hz = abs(max(X(:,3)) - min(X(:,3)))/2;
% 
% delta1 = sqrt(9/(8*hx*hy*hz));
% delta2 = sqrt(27/(8*hx*hy*hz));
% Q = [delta1*q1', delta1*q2', delta1*q3', delta2*q4'];
% P = Q'*Q;
% S = diag(alpha, 0);
% AK = PNstar'*Mtilde*PNstar + (eye(NK)-PN)'*Q*(P\S)/P*Q'*(eye(NK)-PN);