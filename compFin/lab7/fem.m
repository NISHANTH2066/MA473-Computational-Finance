function [ U ] = fem(p,q,f)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    h=0.1;
    x=0:h:1;
    n=length(x);
    A=eye(n);
    F=zeros(n,1);
    for i=2:(n-1)
        A(i,i-1)=-(p(x(i))+p(x(i-1)))/(2*h);
        A(i,i+1)=-(p(x(i))+p(x(i+1)))/(2*h);
        A(i,i)=(2*p(x(i))+p(x(i-1))+p(x(i+1)))/(2*h)+h*q(x(i));
        F(i)=h*f(x(i));
    end
    U=A\F;
    fig=figure();
    plot(x,U);
    hold on;
    
    for i=2:(n-1)
        xi=x(i);x2=x(i+1);x1=x(i-1);
        A(i,i)=(p(x1)+4*p((xi+x1)/2)+p(xi))/(6*h);
        A(i,i)=A(i,i)+(p(x2)+4*p((xi+x2)/2)+p(xi))/(6*h);
        A(i,i)=A(i,i)+h*(q(xi)+q((x1+xi)/2)+q(xi)+q((x2+xi)/2))/6;
        A(i,i-1)=-(p(xi)+4*p((xi+x1)/2)+p(x1))/(6*h);
        A(i,i-1)=A(i,i-1)+h*q((x1+xi)/2)/6;
        A(i,i+1)=-(p(xi)+4*p((xi+x2)/2)+p(x2))/(6*h);
        A(i,i+1)=A(i,i+1)+h*q((xi+x2)/2)/6;
        F(i)=h*(2*f(xi)+2*f((xi+x1)/2)+2*f((xi+x2)/2))/6;
    end
    U=A\F;
    plot(x,U);
end

