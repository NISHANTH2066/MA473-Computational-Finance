function [ x ] = mypsor( A,b)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
    Nmax=500;
    xprev=rand(size(b));
    L=-tril(A,-1);
    U=-triu(A,+1);
    D=A+L+U;
    B=diag(1./diag(A))*(L+U);
    w=2/(1+sqrt(1-norm(B)^2));
    r=b;x=xprev;y=xprev;
    for k=1:Nmax
        for i=1:length(b)
            r(i)=b(i);
            for j=1:i-1
                r(i)=r(i)-A(i,j)*x(j);
            end
            for j=i:length(b)
                r(i)=r(i)-A(i,j)*xprev(j);
            end
            x(i)=max(0,xprev(i)+w*r(i)/A(i,i));
            %y(i)=-r(i)+A(i,i)*(x(i)-xprev(i));
        end
        if norm(x-xprev)<10^(-5)
            break
        end
        xprev=x;
    end
end

