function [ aerr ] = eulerg( Xo,dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    t=0:dt:4;
    n=length(t);
    X=t;S=t;
    X(1)=Xo;
    S(1)=Xo;
    Y=X;
    aerr=0;
    berr=0;
    mu=10;sig=1;
    for it=1:100
        W=randn(2*n,1);
        for i=2:n
            X(i)=X(i-1)-mu*X(i-1)*dt+sig*W(2*i)*sqrt(dt);
            Y(i)=Y(i-1)-mu*Y(i-1)*dt+sig*W(i)*sqrt(dt/2);
        end
        for i=n+1:2*n
            Y(i)=Y(i-1)-mu*Y(i-1)*dt+sig*W(i)*sqrt(dt/2);
        end
        aerr=aerr+abs(X(n)-Y(2*n));
    end
    aerr=aerr/100;
end



