function [ aerr,berr ] = eulermar( Xo,dt )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    t=0:dt:1;
    n=length(t);
    X=t;S=t;
    X(1)=Xo;
    S(1)=Xo;
    Y=X;
    aerr=0;
    berr=0;
    mu=0.75;sig=0.3;
    for it=1:100
        W=sqrt(dt)*randn(n,1);
        W1=sqrt(dt)*randn(n,1);
        for i=2:n
            X(i)=X(i-1)+mu*X(i-1)*dt+sig*X(i-1)*W(i);
            S(i)=S(i-1)*exp((mu-sig^2/2)*dt+sig*W(i));
            Y(i)=Y(i-1)+mu*Y(i-1)*dt+sig*Y(i-1)*W(i)+0.5*sig*Y(i-1)*sig*(W(i)^2-dt);
        end
        aerr=aerr+abs(S(n)-X(n));
        berr=berr+abs(S(n)-Y(n));
    end
    %plot(t,S);
    %hold on;
    %plot(t,X);
    aerr=aerr/100;
    berr=berr/100;
end

