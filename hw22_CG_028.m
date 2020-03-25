function hw22_CG_028

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

