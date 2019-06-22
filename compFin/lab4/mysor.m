function [ x ] = mysor( A,b )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    L=-tril(A,-1);
    U=-triu(A,1);
    D=diag(diag(A));
    Bj=diag(1./diag(A))*(L+U);
    w=2/(1+sqrt(1-norm(Bj)^2));
    Nmax=500;
    xprev=rand(size(b));
    for i=1:Nmax
        x=((1/w-1)*D+U)*xprev+b;
        x=(D/w-L)\x;
        if norm(x-xprev)<10^(-5)
            break
        end
        xprev=x;
    end
end

