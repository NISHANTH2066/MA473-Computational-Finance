% N=10:1000;
% errm=N;erry=N;
% for i=1:length(N)
%     [erry(i),errm(i)]=eulermar(307,1/N(i));
% end
% loglog(1./N,erry);
% hold on;
% loglog(1./N,errm);
% legend('Euler','Milstein');
% 
% eulermarg(307,0.01);

N=100:5:1000;
erry=N;
for i=1:length(N)
    [erry(i)]=eulerg(0,4/N(i));
end
fig=figure();
loglog(1./N,erry);
legend('Error');

egrap(0,1/100);