% Target function
f=@(x) x.^5+x.^3+x;

% Convergence rates
fprintf('Midpoint rule\n');
demo(@midpoint,f,0,6)
fprintf('Trapezoidal rule\n');
demo(@trapezoidal,f,0,6)
fprintf('Simpson rule\n');
demo(@simpson,f,0,6)
fprintf('Monte Carlo rule\n');
demo(@MonteCarlo,f,0,6)
fprintf('Importance sampling rule\n');
demo(@MonteCarloI,f,0,6)

% Comparing trapezoidal rule and Monte Carlo in High dimensional Integral
a=[0,0,0,0,0];
b=[1,1,1,1,1];
g=@(x1,x2,x3,x4,x5) x1.^2+x2.^2+x3.^2+x4.^2+x5.^2;
fprintf(' n results error\n');
fprintf('---------------------------------------------\n');
for s=5:15
    I=0;
    h=(b-a)/s;
    delta=cumprod(h);
    delta=delta(length(a));
    for i=1:s
        for j=1:s
            for k=1:s
                for l=1:s
                    for m=1:s
                        % Define the function
                        I=I+g(a(1)+(i-1)*h(1),a(2)+(j-1)*h(2),a(3)+(k-1)*h(3),a(4)+(l-1)*h(4),a(5)+(m- 1)*h(5))*delta/2;
                        I=I+g(a(1)+i*h(1),a(2)+j*h(2),a(3)+k*h(3),a(4)+l*h(4),a(5)+m*h(5))*delta/2;
                    end
                end
            end
        end
    end
    fprintf('%7d %12.8f %12.8f\n',s^5,I,abs(I-5/3));
end

fprintf(' n results error\n');
fprintf('---------------------------------------------\n');
for s=5:15
    n=s^5; I=0;
    delta=cumprod(b-a); delta=delta(length(a));
    for i=1:n
        p=rand(1,5).*(b-a)+a;
        I=I+g(p(1),p(2),p(3),p(4),p(5));
    end
    I=I/n*delta;
    fprintf('%7d %12.8f %12.8f\n',n,I,abs(I-5/3));
end



% Midpoint rule
function integral= midpoint( f,a,b,n )
    integral=0; h=(b-a)/n;
    for i=1:2:(n-1)
        integral=integral+2*h*f(a+h*i);
    end
end

% Trapezoidal rule
function integral = trapezoidal( f,a,b,n )
    h=(b-a)/n; integral=0;
    for i=1:n
        integral=integral+h*(f(a+(i-1)*h)+f(a+i*h))/2;
    end
end

% Simpson?s rule
function integral= simpson( f,a,b,n )
    h=(b-a)/n; sum=f(a)+f(b);
    for i=2:2:n-2
        sum=sum+2*f(a+i*h);
    end
    for i=1:2:n-1
        sum=sum+4*f(a+i*h);
    end
    integral=h/3*sum;
end

% Monte Carlo
function I = MonteCarlo( f,a,b,n )
    p=rand(n,1)*(b-a)+a;
    I=sum(f(p)*(b-a)/n);
end

% Monte Carlo Importance sampling
function I = MonteCarloI( f,a,b,n )
    p = sqrt(rand(n,1)*(b^2-a^2)+a^2); % p ~ p(x?= 2*x/(b^2-a^2) distribution
    I = f(p).*(b^2-a^2)./p/2; % f(x)/p(x)
    I = mean(I);
end

% Demo function
function [ ] = demo( method,f,a,b )
    truevalue=integral(f,a,b); %integral is a built-in function of matlab. %This gives the true value of the integral being estimated.
    fprintf(' n h results error\n'); fprintf('------------------------------------------------------\n');
    for i=1:7
        n=10^i;
        value=method(f,a,b,n); %We select a method
        fprintf('%7d %12.7f %12.8f %12.8f\n',n,(b-a)/n,value,abs(value-truevalue));
    end
end
        
        
        
