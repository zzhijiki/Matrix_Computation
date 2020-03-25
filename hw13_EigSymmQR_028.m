
%
% 你的学号和姓名
% 
% 带 Wilkinson 位移的 QR 算法
%
% 请将文件名和函数名中的 xxx 改为你的学号最后三位
%
function hw13_EigSymmQR_028

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
if n==1;
    Eig=Diag;
    return  
    
end
for iter = 1:iter_max
    
    shit = Eig2by2(A(n-1,n),A(n-1,n-1),A(n,n));
    [b,p] = min ( abs( shit-A(n,n) )  ) ;%第一个 变量没用
    shift=shit(p); 
    [c,s] = Givens(A(1,1)-shift , A(2,1));
    G=eye(n);
    G(1:2,1:2)=[[c,-s];[s,c]];
    A=G'*A*G;
    
        
    for i = 1:n-2
        G=eye(n);
        [c,s] = Givens(A(i+1,i), A(i+2,i));
        G(i+1:i+2,i+1:i+2)=[[c,-s];[s,c]];
        A=G'*A*G    ;
    end
    
    Subdiag=diag(A,1);
    index=find(abs(Subdiag)<tol);
    if length(index)>0
       break
    end
    
end
p=index(1);

C=A(1:p,1:p);
D=A(p+1:n,p+1:n);


Eig=[ EigSymmQR( diag(C),diag(C,1) );  EigSymmQR( diag(D),diag(D,1) )  ];%递归算法





%%%%% 子函数：Givens 变换
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

%%%%% 子函数：Eig2by2 -- 计算 2 阶对称矩阵的特征值
function Eig = Eig2by2(a,b,c)
%
% A = [b, a; a, c]
%
Sqrt_Delta = sqrt((b-c)*(b-c)+4*a*a);
Eig = [(b+c-Sqrt_Delta)/2; (b+c+Sqrt_Delta)/2];
