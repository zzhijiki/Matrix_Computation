% 
% GMRES
%

function hw22_GMRES_028

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
