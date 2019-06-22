function [ U ] = thetam(e)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
    tim=1;K=10;
    r=0.06;sig=0.3;del=0;
    dx=1;
    dt=dx^2/20;
    Smax=20;
    S=0:dx:Smax;
    T=0:dt:tim;
    m=length(S);
    n=length(T);
    U=zeros(m,n);
    A=zeros(m);
    B=zeros(m);
    U(:,n)=max(S-K,0);
    U(1,:)=0;
    U(m,:)=Smax-K*exp(-r*(tim-T));
    aa = @(x) 0.5*sig^2*x^2/dx^2;
    bb = @(x) 0.5*(r-del)*x/dx;
    for j=n-1:-1:1
        A(1,1)=1;A(m,m)=1;B=A;
        for i=2:m-1
            A(i,i+1)=(aa(S(i))+bb(S(i)))*e;
            A(i,i-1)=(aa(S(i))-bb(S(i)))*e;
            A(i,i)=e*(-2*aa(S(i))-r)+1/dt;
            B(i,i+1)=-(1-e)*(aa(S(i))+bb(S(i)));
            B(i,i-1)=-(1-e)*(aa(S(i))-bb(S(i)));
            B(i,i)=-(1-e)*(-2*aa(S(i))-r)+1/dt;
        end
        temp=U(:,j+1);
        U(:,j)=B\(A*temp);
        U(1,:)=0;
        U(m,:)=Smax-K*exp(-r*(tim-T));
    end
    fig=figure();
    surf(S,T,U');
    fig=figure();
    plot(S,U(:,1));
    hold on;
    plot(S,U(:,end));
end

