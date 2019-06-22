function [ x ] = mygauss( A,b )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
    L=-tril(A,-1);
    U=-triu(A,1);
    D=diag(diag(A));
    xprev=rand(size(b));
    Nmax=500;
    for i=1:Nmax
        x=U*xprev+b;
        x=(D-L)\x;
        if norm(x-xprev)<10^(-5)
            break
        end
        xprev=x;
    end
end

