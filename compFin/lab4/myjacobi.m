function [ x ] = myjacobi(A,b)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    L = -tril(A,-1);
    U = -triu(A,1);
    D = diag(diag(A));
    xprev=rand(size(b));
    for i=1:500
        x=(L+U)*xprev+b;
        x=D\x;
        if norm(x-xprev)<10^(-5)
            break
        end
        xprev=x;
    end
end

