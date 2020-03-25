# 矩阵计算课后题MATLAB代码作业

### 第二章

#### 练习2.15

解方程，对应p43的思考题，pass

#### 练习2.16

trick：事实上，追赶法就是LU分解，但是他不需要生成（还原）原来的矩阵，只需要三条对角线就可以了。

对应算法2.11（p55），只需要初始化一些东西

```matlab
练习2.16 代码

function LU_tridiag

%使用自己做的例子
%A=[1,2,0;1,4,6;0,7,9]
%f=[2,3,4]'
%x=[3/4,5/8,-1/24]'    #我们的目标

a=[1,7];
b=[1,4,9];
c=[2,6];
f=[2,3,4];

n=length(b);
alpha=zeros(n,1);  %初始化一些东西
beta=zeros(n-1,1);
y=zeros(n,1);
x=zeros(n,1);        %初始化一些东西


alpha(1)=b(1);     %算法2.11照抄
beta(1)=c(1)/b(1);
y(1)=f(1)/b(1);
for i=2:n-1
    alpha(i)=b(i)-alpha(i-1)*beta(i-1);
    beta(i)=c(i)/alpha(i);
    y(i)=(f(i)-a(i-1)*y(i-1))/alpha(i);
end
alpha(n)=b(n)-a(n-1)*beta(n-1);
y(n)=(f(n)-a(n-1)*y(n-1))/alpha(n);
x(n)=y(n);
for i =n-1:-1:1
    x(i)=y(i)-beta(i)*x(i+1);
end
x
```

#### 练习2.17

感觉就是生成带状矩阵（给你几条对角线），然后做LU分解。

重点就是生成带状矩阵，用个for循环即可

```matlab
练习 2.17代码

n=6; %生成带状矩阵的方法，我们有的东西是 几个对角线 或 次对角线
bl=2;
bu=3;
A=diag(diag(rand(n,n)));
for k=1:bl
    A=A+diag(rand(n-k,1),-k);
end
for k=1:bu
    A=A+diag(rand(n-k,1),k);
end
A  
```

### 第三章

```matlab
练习3.29代码:

R=[1,2,3,4,5;0,5,6,7,8;0,0,8,9,10;0,0,0,10,11;0,0,0,0,4];  %构造了一个题目
v=[1,2,3,4,5]';
u=[4,5,6,7,8]';
u*v';
A=R+u*v'   %构造了一个题目

n=length(v);
P=eye(n);  %记录用过的正交矩阵
for i = n:-1:3
    G=eye(n);
    G1=givens(A(i-1,1),A(i,1));  %matlab 内置函数givens(),参数两个数，返回2*2矩阵
    G(i-1:i,i-1:i)=G1;
    A=G*A;                      %做givens变换
    P=P*G';                   %记录用过的正交矩阵
end    
A   %验证是否为上HESS
P*A %验证回到了原来的矩阵


Q=eye(n);   %对上HESS做QR
for k=1:n
    for i=k+1:n
    G=eye(n);
    G1 = givens(A(k,k), A(i,k));   %matlab 内置函数givens(),参数两个数，返回2*2矩阵

    B=[A(k,k:n);A(i,k:n)];   %做givens变换
    B=G1*B;                     %做givens变换
    A(k,k:n)=B(1,:);          %做givens变换
    A(i,k:n)=B(2,:);         %做givens变换
    
    C=[Q(1:n,k),Q(1:n,i)];    %记录用过的正交矩阵
    C=C*G1';                 %记录用过的正交矩阵
    Q(1:n,k)=C(:,1);        %记录用过的正交矩阵
    Q(1:n,i)=C(:,2);         %记录用过的正交矩阵
    
    end
end
A  %验证是否为上三角
P*Q*A %验证回到了原来的矩阵
```

#### 练习3.30

列主员QR，主要针对非满秩情况，因为他觉得，前面一串不好看，所以要换一下列。

思路: 基于Householder的QR分解  + 列交换（采用相关代码的程序PLU.m）

#### 难点：

vecnorm函数：求出矩阵每一列的2范数，列交换是基于范数大小交换的

代码3.30 ：将两个代码有机组合，成功了（R是上三角的，e-16归零的话。 且前2（=rank）列有东西，后一列为0），于是没检查，可能不对。

```matlab
A=[0,2,4;0,3,2;0,1,5;0,6,7];  %矩阵A肯定要非满秩

p = 1 : n; 
[m,n]=size(A);
Q=eye(m);
for k = 1:n
    B=A(k:n,k:n);
    [a_max,l] = max(abs(vecnorm( B )));   %求出矩阵每一列的2范数

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
Q*R  %验证对不对

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
```

#### 练习3.31

故弄玄虚，其实是对A转置求QR分解

LQ分解：

[Q1,R1] = qr (A');

A=(Q1 *R1)'=R1'* Q1'=L* Q;

所以下三角阵L=R1',正交阵Q=Q1'

### 第四章

#### 练习 4.15

对应算法4.7 ，需要调用相关代码的householder函数

```matlab
A=[1,2,3,4;5,6,7,8;2,3,4,5;3,4,5,6]
hess(A)     %与matlab内置函数比较
[n,n]=size(A);
Q=eye(n);

for k=1:n-2    %抄书
    [beta,v]=House(A(k+1:n,k)) ;
    A(k+1:n,k:n)=A(k+1:n,k:n)-beta*v*(v'*A(k+1:n,k:n));
    A(1:n,k+1:n)=A(1:n,k+1:n)-beta*A(1:n,k+1:n)*v*v';
    Q(k+1:n,:)=Q(k+1:n,:)-beta*v*(v'*Q(k+1:n,:));
end
A  %这是我们抄书得到的东西，好像差了e-16次的误差，应该没关系（若要修改的话，加一个tol，归0即可）
Q   %需要返回的Q
Q'*A*Q  %发现回到了A，验证算法正确性


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
```

### 第五章

#### 练习5.10

和上面练习 4.15的HESS没有本质区别

#### 练习5.11

对应homework

递归算法

练习5.10代码

In [ ]:

```matlab
function hw13_EigSymmQR_0学号

global iter_max tol
iter_max = 10; % 计算每个特征值的最大迭代步数
tol = 1e-12;   % 相对精度

% 生成 n 阶随机对称三对角矩阵
n = 100;
Diag = randn(n,1);      % 主对角线
SubDiag = randn(n-1,1); % 次对角线

% 计算特征值
Eig = EigSymmQR(Diag,SubDiag);

% 测试
A = diag(Diag) + diag(SubDiag,1) + diag(SubDiag,-1);

TrueEig = eig(A);  % 用 MATLAB 自带函数计算特征值
fprintf('relerr = %.4e\n', norm(sort(Eig)-sort(TrueEig)) / norm(TrueEig))


%%%%% 带 Wilkinson 位移的对称 QR 算法
function Eig = EigSymmQR(Diag,SubDiag)
%
% Input:  Diag -- diagonal vector of the symmetric tridiagonal matrix
%         SubDiag -- subdiagonal vector
% Output: Eig -- eigenvalues
%
global iter_max tol

%
% 你的代码
%

A = diag(Diag) + diag(SubDiag,1) + diag(SubDiag,-1);
n=length(Diag);
if n==1;   %终止条件
    Eig=Diag;
    return  
    
end
for iter = 1:iter_max
    
    shit = Eig2by2(A(n-1,n),A(n-1,n-1),A(n,n));
    [b,p] = min ( abs( shit-A(n,n) )  ) ;%第一个 变量没用
    shift=shit(p); 
    [c,s] = Givens(A(1,1)-shift , A(2,1));    %given
    G=eye(n);   %given
    G(1:2,1:2)=[[c,-s];[s,c]];     %given
    A=G'*A*G;  %第一次QR迭代
    
        
    for i = 1:n-2   %QR迭代
        G=eye(n);
        [c,s] = Givens(A(i+1,i), A(i+2,i));
        G(i+1:i+2,i+1:i+2)=[[c,-s];[s,c]];
        A=G'*A*G    ;
    end       %QR迭代
    
    Subdiag=diag(A,1);
    index=find(abs(Subdiag)<tol);    %找到我的分割点
    if length(index)>0
       break
    end
    
end
p=index(1);   %找到我的分割点中最前面的那个

C=A(1:p,1:p);
D=A(p+1:n,p+1:n);


Eig=[ EigSymmQR( diag(C),diag(C,1) );  EigSymmQR( diag(D),diag(D,1) )  ];%递归算法





%%%%% 子函数:Givens 变换
function [c, s] = Givens(a, b)
%
% Input:  a -- REAL scalar
%         b -- REAL scalar
% Output: c -- REAL rotation matrix element
%         s -- REAL rotation matrix element
%
if ( b == 0.0 )
    c = 1.0; s = 0.0;
elseif ( abs(b) > abs(a) )
    tau = a / b;
    s = 1.0 / sqrt( 1.0 + tau^2 );
    c = tau * s;
else
    tau = b / a;
    c = 1.0 / sqrt( 1.0 + tau^2 );
    s = tau * c;
end

%%%%% 子函数:Eig2by2 -- 计算 2 阶对称矩阵的特征值
function Eig = Eig2by2(a,b,c)
%
% A = [b, a; a, c]
%
Sqrt_Delta = sqrt((b-c)*(b-c)+4*a*a);
Eig = [(b+c-Sqrt_Delta)/2; (b+c+Sqrt_Delta)/2];
```

### 第六章

#### 练习6.16

possion 方程的红黑排序 SOR 和SSOR 迭代

书上给出的算法6.4 是红黑排血的G-S迭代

所以我们要写出 SOR 的迭代格式（算法6.5只是一维的，而我们红黑排序必然是二维的，所以要做修改）

即

$SOR:u_{i,j}=(1−w)∗u_{i,j}+(w/4)∗[u_{i-1,j}+u_{i,j-1}+u_{i,j}+1+h∗h∗f]$

这就是没有进行红黑排序的SOR，其实我们有这个算法的源代码Poisson_SOR_omega，所以只要在这个基础上修改成红黑排序即可

而SSOR与SOR公式是一样的，他只是多了个循环（算法6.6）

或者可以参照相关代码Poisson_SSOR_omega和Poisson_SOR_omega 可以看出区别

```matlab
所谓红黑排序只要找到红色的点和黑色的点就可以了
首先在老师给的代码中的循环有
i=2:N+1
j=2:N+1

那么, 
红色点特征:i+j=2p(偶数)=[4,6...,2N+2]
黑色点特征:i+j=2p+1(奇数)=[5,7...,2N+1]

所以关于红色,只要在构造一个p做循环就可以了
for i = 2:N+1
    for p = 2:N+1
        j=2*p-i;
        
可以看到j+i =2p,符合要求,但是j的取值要2~N+1,加一个判断条件即可  if j<N+2 && j>1 

所以关于黑色,只要在构造一个p做循环就可以了,这里p只要到N就可以了 
for i = 2:N+1
    for p = 2:N
        j=2*p+1-i;
       
        
可以看到j+i =2p+1,符合要求,但是j的取值要2~N+1,加一个判断条件即可  if j<N+2 && j>1 
```

代码练习6.16 ，参考Poisson_SOR_omega ，除了循环一个字也没改

```matlab
% 测试 SOR 对参数 omega 的敏感性
%
%

clear all;
close all;
N = 8;
h = 1/(N+1); % step size

% initial guess
v0 = zeros(N+2,N+2);

% boundary condition
x = 0:h:1;  y = 0:h:1;
v0(1,:) = 0.25 * (y .* y);
v0(N+2,:) = 0.25 * (1 + y .* y);
v0(:,1) = 0.25 * (x .* x);
v0(:,N+2) = 0.25 * (x .* x + 1);

f = -1; % right-hand side f(x,y)

% true solution
[X,Y] = meshgrid(x);
vt = 0.25*(X.*X + Y.*Y);
norm_vt = norm(vt(:));

tol = 1e-6; % stopping criteria
iter_max = 1000; % maximun number of iterations

hh = h*h;  % h^2

% values of omega
Omega = 1.2:0.01:1.7;

Iter = zeros(size(Omega));  % iteration number for different omega

% start SOR
for k_omega = 1 : length(Omega)
    omega = Omega(k_omega);
    v = v0;
    for k = 1 : iter_max
        for i = 2:N+1
            for p = 2:N+1
                j=2*p-i;
                
                if j<N+2 && j>1 
                    vij_old = v(i,j);
                    v(i,j) = 0.25* (hh*f + v(i+1,j) + v(i-1,j)+ v(i,j+1) + v(i,j-1));
                    v(i,j) = (1-omega)*vij_old + omega*v(i,j);
                end
            end
        end
        
        for i = 2:N+1
            for p = 2:N
                j=2*p+1-i;
                if j<N+2 && j>1 
                    vij_old = v(i,j);
                    v(i,j) = 0.25* (hh*f + v(i+1,j) + v(i-1,j) ...
                    + v(i,j+1) + v(i,j-1));
                    v(i,j) = (1-omega)*vij_old + omega*v(i,j);
                end    
            end
        end
        relerr = norm(v(:) - vt(:)) / norm_vt;
        if relerr < tol
            break;
        end
    end
    Iter(k_omega) = k;
end

% plot the result
plot(Omega,Iter,'o-');
axis([1.2,1.7,20,70])
xlabel('\omega'); ylabel('iteration number');
legend([int2str(N),' x ', int2str(N)]);
title('SOR with different \omega');
```

### 第七章

#### 练习 7.4 GMRES方法

对应课后作业，以及算法7.5照抄（p291）

难点：

1. 计算离散的方程，二维，一维possion公式，有了kron 系数矩阵容易构造
2. 求导后的右端项所有乘法都是 .* ，并且还要meshgrid计算normb之类的乱七八糟的东西，记在书上。
3. 傻瓜作图，semilogy(relres);

```matlab
function hw22_GMRES_0学号

global N hx hy

% initialization
N = 64; n = N*N;
a = 0; b = 1;
hx = (b-a)/(N+1);
hy = hx;
tol = 1e-6;
IterMax = 200;
u0 = zeros(N*N,1);   % initial guess

% 右端项
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%
x=0:hx:1; x=x(2:N+1);
[X,Y]=meshgrid(x);
f=(3-2.*X).*(1-Y).*Y+(3-2.*Y).*(1-X).*X+(1-X).*X.*(1-Y).*Y;  f=f(:);
b=hx*hx*f;
norm_b=norm(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize workspace
e1 = zeros(n,1);  e1(1) = 1.0;
cs = zeros(IterMax,1);  
sn = zeros(IterMax,1);

% begin GMRES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%
V=b/norm_b;
e=norm_b*e1;
for j=1:IterMax
    w=Ax(V(:,j));
    for i=1:j
        H(i,j)=V(:,i)'*w;
        w=w-H(i,j)*V(:,i);
    end
    H(j+1,j)=norm(w);
    if H(j+1,j)==0
        m=j; break
    end
    V(:,j+1)=w/H(j+1,j);
    for i=1:j-1  
        H([i,i+1],j)=[cs(i) sn(i); -sn(i) cs(i)]*H([i,i+1],j);
    end
    [cs(j),sn(j)]=mygivens(H(j,j),H(j+1,j));
    H(j,j)=cs(j)*H(j,j)+sn(j)*H(j+1,j);
    H(j+1,j)=0;
    e(j:j+1)=[cs(j) sn(j); -sn(j) cs(j)]*[e(j),0]';
    relres(j)=norm(e(j+1))/norm_b;   %QR分解算res
    if relres(j)<tol
        m=j;break
    end    
end
m=j;
y=H(1:m,1:m)\e(1:m);
u=u0+V(:,1:m)*y;
fprintf('Iter=%d, relres=%.4e\n',m, relres(m));
semilogy(relres);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Givens 变换
function [c, s] = mygivens(a, b)
% compute the Givens rotation matrix parameters for a and b.
%
% Input
%   a -- REAL scalar
%   b -- REAL scalar
% Output
%   c -- REAL rotation matrix element
%   s -- REAL rotation matrix element

if ( b == 0.0 )
    c = 1.0;
    s = 0.0;
elseif ( abs(b) > abs(a) )
    tau = a / b;
    s = 1.0 / sqrt( 1.0 + tau^2 );
    c = tau * s;
else
    tau = b / a;
    c = 1.0 / sqrt( 1.0 + tau^2 );
    s = tau * c;
end

% 计算系数矩阵与向量的乘积
function y = Ax(u)
global N hx hy
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%
v1=ones(1,N);   %  Tn 的 diag
v2=ones(1,N-1);  %  Tn 的 次 diag
T=diag(2*v1)+diag(-1*v2,1)+diag(-1*v2,-1);   %  Tn
I=eye(N);  % I
D=diag(0.5*v2,1)+diag(-0.5*v2,-1);  %  D

A=kron(I,T)+kron(T,I)+hx.*kron(I,D)+hx.*kron(D,I)+(hx*hx).*kron(I,I);  %A

y=A*u;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ===== End =====
```

#### 练习7.5 CG算法：

对应227页 算法7.7（套用上一题的模板，把GMRES改成了CG）

```matlab
代码7.5

function hw22_CG_0学号

global N hx hy A

v1=ones(1,N);   %  Tn 的 diag
v2=ones(1,N-1);  %  Tn 的 次 diag
T=diag(2*v1)+diag(-1*v2,1)+diag(-1*v2,-1);   %  Tn
I=eye(N);  % I
D=diag(0.5*v2,1)+diag(-0.5*v2,-1);  %  D

A=kron(I,T)+kron(T,I)+hx.*kron(I,D)+hx.*kron(D,I)+(hx*hx).*kron(I,I);  %A
% initialization
N = 64; n = N*N;
a = 0; b = 1;
hx = (b-a)/(N+1);
hy = hx;
tol = 1e-6;
IterMax = 1500;
u0 = zeros(N*N,1);   % initial guess

% 右端项
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%
x=0:hx:1; x=x(2:N+1);
[X,Y]=meshgrid(x);
f=(3-2.*X).*(1-Y).*Y+(3-2.*Y).*(1-X).*X+(1-X).*X.*(1-Y).*Y;  f=f(:);
b=hx*hx*f;
norm_b=norm(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize workspace
e1 = zeros(n,1);  e1(1) = 1.0;
cs = zeros(IterMax,1);  
sn = zeros(IterMax,1);

% begin GMRES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%
r=b;
beta=norm_b;
roll0=r;
for j=1:IterMax
    roll=r'*r;
    if j>1
        u=roll/roll0;
        p=r+u*p;
    else
        p=r;
    end
    q=A*p;
    epsilon=roll/(p'*q);
    x=x+epsilon*p;
    r=r-epsilon*q;
    relres(j)=norm(r)/beta;

    if relres(j)<tol
        
        break
    end
    roll0=roll;
end

fprintf('Iter=%d, relres=%.4e\n',j, relres(j));
semilogy(relres);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Givens 变换
function [c, s] = mygivens(a, b)
% compute the Givens rotation matrix parameters for a and b.
%
% Input
%   a -- REAL scalar
%   b -- REAL scalar
% Output
%   c -- REAL rotation matrix element
%   s -- REAL rotation matrix element

if ( b == 0.0 )
    c = 1.0;
    s = 0.0;
elseif ( abs(b) > abs(a) )
    tau = a / b;
    s = 1.0 / sqrt( 1.0 + tau^2 );
    c = tau * s;
else
    tau = b / a;
    c = 1.0 / sqrt( 1.0 + tau^2 );
    s = tau * c;
end
```

