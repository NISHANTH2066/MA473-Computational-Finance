function [ U ] = theta(a)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    % a is theta
    dx=0.1;dt=0.001;
    xmax=1;
    xmin=-5;
    r=0.06;sig=0.3;del=0;K=10;
    tau=sig^2/2;
    dx=0.05;
    dt=dx^2/2;
    k=dt/dx/dx;
    X=xmin:dx:xmax;
    T=0:dt:tau;
    m=length(X);
    n=length(T);
    A=zeros(m);
    B=zeros(m);
    q=2*r/sig^2;
    qd=2*(r-del)/sig^2;
    
    for i=2:m-1
        A(i,i-1)=k*a;
        B(i,i-1)=-k*(1-a);
        A(i,i)=-2*k*a-1;
        B(i,i)=2*k*(1-a)-1;
        A(i,i+1)=k*a;
        B(i,i+1)=-k*(1-a);
    end
    A(1,1)=1;
    B(1,1)=1;
    B(m,m)=1;
    A(m,m)=1;
    U=zeros(m,n);
    % U(1,:) is already zero
    U(:,1)=init(X,qd);
    U(m,:)=bc2(xmax,qd,T);
    U(:,1)=init(X,qd);
    for j=2:n
        temp=U(:,j-1);
        temp(1)=U(1,j);
        temp(end)=U(end,j);
        U(:,j)=A\(B*temp);
    end
    fig=figure();
    surf(U);
    
    V=zeros(m,n);
    for i=1:m
        for j=1:n
            aa=-(qd-1)*X(i)/2;
            bb=-(((qd-1)^2)/4+q)*T(j);
            V(i,j)=K*exp(aa+bb)*U(i,j);
        end
    end
    fig=figure();
    S = K*exp(X);
    tt = 1-T/sig;
    surf(X,T,flip(V'));
    
    fig=figure();
    plot(S,V(:,end));
    hold on;
    plot(S,V(:,1));
end

function [ X ] = init(x,qd)
    X = exp(x.*(qd+1)/2)-exp(x.*(qd-1)/2);
    X = max(X,0);
end

function [ X ] = bc2(x,qd,t)
    X = exp(x*(qd+1)/2+t.*(qd+1)^2/4);
end

