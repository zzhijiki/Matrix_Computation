A=[1,2,3,4;5,8,7,8;2,3,4,5;3,4,7,6]
hess(A)
[n,n]=size(A);
Q=eye(n);

for k=1:n-2
    [beta,v]=House(A(k+1:n,k)) ;
    A(k+1:n,k:n)=A(k+1:n,k:n)-beta*v*(v'*A(k+1:n,k:n));
    A(1:n,k+1:n)=A(1:n,k+1:n)-beta*A(1:n,k+1:n)*v*v';
    Q(k+1:n,:)=Q(k+1:n,:)-beta*v*(v'*Q(k+1:n,:));
end
A
Q'*A*Q


function [beta,v]=House(x)  %调用网站上给的householder函数
n = length(x);
sigma = dot(x(2:n),x(2:n));
v = x;
if sigma == 0 
    if x(1) < 0
        v(1) = 2*x(1);
        beta = 2/(v(1)*v(1));
    else
        v(1) = 1; 
        beta = 0;
    end
else
    alpha = sqrt(x(1)*x(1)+sigma);
    if x(1) < 0
        v(1) = x(1) - alpha;
    else
        v(1) = -sigma/(x(1)+alpha);
    end
    beta = 2/(v(1)*v(1)+sigma);
end
end