function [ V ] = thetashi(e)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Tmax=1;r=0.04;sig=0.25;del=0.1;K=10;
    dx=0.1;dt=0.01;
    T=0:dt:Tmax;T=Tmax-T;
    Smax=20; % Not needed
    X=0:dx:1;
    m=length(X);
    n=length(T);
    A=eye(m);
    B=eye(m);
    a = @(x) 0.5*sig^2*(x*(1-x))^2/dx^2;
    b = @(x) 0.5*(r-del)*x*(1-x)/dx;
    c = @(x) r*(1-x)+del*x;
    for i=2:m-1
        A(i,i)=e*(-2*a(X(i))-c(X(i)))-1/dt;
        A(i,i-1)=e*(a(X(i))-b(X(i)));
        A(i,i+1)=e*(a(X(i))+b(X(i)));
        B(i,i-1)=-(1-e)*(a(X(i))-b(X(i)));
        B(i,i+1)=-(1-e)*(a(X(i))+b(X(i)));
        B(i,i)=(1-e)*(2*a(X(i))+c(X(i)))-1/dt;
    end
    U=zeros(m,n);
    % Already done U(1,:)
    U(end,:)=exp(-del*T);
    U(:,1)=max(2*X-1,0);
    for j=2:n
        temp=U(:,j-1);
        %U(:,j)=myjacobi(A,B*temp);
        %U(:,j)=A\(B*temp);
        %U(:,j)=mygauss(A,B*temp);
        U(:,j)=mysor(A,B*temp);
        U(end,:)=exp(-del*T);
    end
    V=zeros(m,n);
    for i=1:m
        for j=1:n
            V(i,j)=K*U(i,j)/(1-X(i));
        end
    end
    fig=figure();
    surf(X*K./(1-X),T,V'); % Note that V(m,:)=Inf
    fig=figure();
    plot(X*K./(1-X),V(:,end));
    hold on;
    plot(X*K./(1-X),V(:,1));
end

