function [ U ] = atheta( e )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    T=1;K=10;r=0.25;sig=0.6;del=0.2;
    qd=2*(r-del)/sig^2;q=2*r/sig^2;
    dx=0.1;dt=dx^2/2;
    lam=dt/dx^2;
    tau=T*sig^2/2;
    t=0:dt:tau;
    xmin=-4;xmax=1;
    X=xmin:dx:xmax;
    m=length(X);
    n=length(t);
    C=zeros(m,m);
    for i=2:m-1
        C(i,i)=1+2*lam*e;
        C(i,i-1)=-lam*e;
        C(i,i+1)=-lam*e;
    end
    C=C(2:m-1,2:m-1);
    g = @(x,t) exp(t/4*((qd-1)^2+4*q))*max(exp(x*(qd-1)/2)-exp(x*(qd+1)/2),0);
    U=zeros(m,n);
    U(1,:)=g(xmin,t);
    U(m,:)=g(xmax,t);
    U(:,1)=g(X,0);
    for j=2:n
        b=X(2:end-1)';
        for ii=2:length(b)-1
            i=ii+1;
            b(ii)=U(i,j-1)+lam*(1-e)*(U(i+1,j-1)-2*U(i,j-1)+U(i-1,j-1));
        end
        b(1)=U(2,j-1)+lam*(1-e)*(U(3,j-1)-2*U(2,j-1)+U(1,j-1))+e*lam*U(1,j);
        b(end)=U(m-1,j-1)+lam*(1-e)*(U(m,j-1)-2*U(m-1,j-1)+U(m-2,j-1))+e*lam*U(m,j);
        gg=b;
        for ii=1:length(b)
            i=ii+1;
            gg(ii)=g(X(i),t(j));
        end
        U(2:end-1,j)=mypsor(C,b-C*gg)+gg;
    end
    V=U;
    for i=1:m
        for j=1:n
            V(i,j)=K*exp(-(qd-1)*X(i)/2-((qd-1)^2/4+q)*t(j))*U(i,j);
        end
    end
    fig=figure();
    surf(K*exp(X),(T-t)/tau,V');
    fig=figure();
    plot(K*exp(X),V(:,end));
    hold on;
    plot(K*exp(X),V(:,1));
end

