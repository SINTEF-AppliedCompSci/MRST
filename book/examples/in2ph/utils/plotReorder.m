% Script to reorder and plot the diagonal blocks of a Jacobi matrix. This
% can e.g., be called after line 103 has been executed in the private
% function newtonRaphson2ph from the incomp module

[p,q,r,s]=dmperm(J);
figure;
spy(J(p,q));
n = find(diff(r)>1);
hold on
for j=n
    plot(r([j j+1 j+1 j j])-.5,r([j j j+1 j+1 j])-.5,'-r');
end
hold off