A=[0,2,4;0,3,2;0,1,5;0,6,7];

p = 1 : n; 
[m,n]=size(A);
Q=eye(m);
for k = 1:n
    B=A(k:n,k:n);
    [a_max,l] = max(abs(vecnorm( B )));

    l = l+k-1;
    if l~=k
        t = A(:,k); A(:,k) = A(:,l); A(:,l) = t;
        tmp = p(k); p(k) = p(l); p(l) = tmp;
    end
    x=A(k:m,k);
    [b,v]=House(x);
    A(k:m,k:n)= A(k:m,k:n)-b*v*(v'*A(k:m,k:n));
    Q(:,k:m)=Q(:,k:m)-b*(Q(:,k:m)*v)*v'  ;  
end
R=A
Q
Q*R

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
