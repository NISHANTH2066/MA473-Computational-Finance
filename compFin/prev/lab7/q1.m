function [ U ] = q1()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    dx = 0.1;
    h = dx;
    X = 0:dx:1;
    n = length(X);
    I1 = zeros(n);
    I2 = zeros(n);
    F = zeros(n,1);
    for i=1:n
        if i==1
            
        elseif i==n
            
        else
            I1(i,i)=(p(X(i-1))+2*p(X(i))+p(X (i+1)))/(2*h);
            I1(i,i-1)=-(p(X(i-1))+p(X(i)))/(2*h);
            I1(i,i+1)=-(p(X(i+1))+p(X(i)))/(2*h);
            I2(i,i)=h*q(X(i));
            F(i) = h*f(X(i));
        end
    end
    A = I1+I2;
    A = A(2:n-1,2:n-1);
    %disp(A);
    %disp(F);
    F = F(2:n-1);
    U = A\F;
    U = [0; U ; 0];
    fig = figure();
    plot(X,U);
    grid on;
    %title('Trapezoidal Rule');
    saveas(fig,'q1_trap.jpg');
    
    I1 = zeros(n);
    I2 = zeros(n);
    F = zeros(n,1);
    % Simsons Rule
    for i=1:n
        if i==1
            
        elseif i==n
            
        else
            t1 = (X(i-1)+X(i))/2;
            I1(i,i)=(p(X(i-1))+4*p(t1)+2*p(X(i))+4*p(t1+h)+p(X (i+1)))/(6*h);
            %I1(i,i)=2*(p(X(i-1))+4*p(X(i))+p(X(i+1)))/(6*h);
            I1(i,i-1)=-(p(X(i-1))+4*p(t1)+p(X(i)))/(6*h);
            I1(i,i+1)=-(p(X(i+1))+4*p(t1+h)+p(X(i)))/(6*h);
            %I2(i,i)=h*4*q(X(i))/3;
            I2(i,i)=h*(q(t1)+2*q(X(i))+q(t1+h))/6;
            %F(i) = h*4*f(X(i))/3;
            F(i) = h*(2*f(X(i))+2*f(t1)+2*f(t1+h))/6;
        end
    end
    A = I1+I2;
    A = A(2:n-1,2:n-1);
    F = F(2:n-1);
    U = A\F;
    U = [0; U ; 0];
    hold on;
    plot(X,U);
    legend('trapezoidal','simpson');
%     fig = figure();
%     plot(X,U);
%     grid on;
%     title('Simsons Rule');
%     saveas(fig,'q1_simp.jpg');

    % Or do it by taking A(1,1)=A(N,N)=1
end

function[y] = p(x)
    y=x+1;
end

function[y] = q(x)
    y = x*x+2;
end

function[y] = f(x)
    y = x*x-4;
end
