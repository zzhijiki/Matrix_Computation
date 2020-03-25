%
% 你的学号和姓名
% 
% 用三种方法实现 QR 分解，并分析其稳定性
%
% 请将文件名和函数名中的 xxx 改为你的学号最后三位
%
function hw12_QRStability_028

m = 400; 
n = 300;

for cnd = 10.^[0:6]  % 条件数
    for k = 1 : 10   % 每个条件数测试 10 次
        A = randn(m,n);
        [U,S,V] = svd(A);  % 生成两个随机的酉矩阵
        
        % singular values range from 1 to cnd, with uniformly distributed logarithms
        sd = [1,cnd, exp(rand(1,n-2)*log(cnd))];
        A = U(:,1:n)*diag(sd)*V';   % 生成指定奇异值的随机矩阵
        
        [Q_GS,R] = QR_GS(A);
        err_GS(k) = norm(A-Q_GS*R)/cnd;
        err_GS_orth(k) = norm(Q_GS'*Q_GS - eye(n));
        
        [Q_MGS,R] = QR_MGS(A);
        err_MGS(k) = norm(A-Q_MGS*R)/cnd;
        err_MGS_orth(k) = norm(Q_MGS'*Q_MGS - eye(n));
        
        [Q_H,R] = QR_House(A);
        err_H(k) = norm(A-Q_H*R)/cnd;
        err_H_orth(k) = norm(Q_H'*Q_H - eye(m));
        
    end
    
    %%% 输出测试结果
    fprintf('\ncnd=%d\n',cnd);
    fprintf('   GS: relerr = %.4e, err_orth=%.4e\n', max(err_GS), max(err_GS_orth))
    fprintf('  MGS: relerr = %.4e, err_orth=%.4e\n', max(err_MGS),max(err_MGS_orth))
    fprintf('House: relerr = %.4e, err_orth=%.4e\n', max(err_H),  max(err_H_orth))
end

%%% 用 GS 方法实现 QR 分解
function [Q,R] = QR_GS(A)
%
% 你的代码
%
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n);
R(1,1)=norm(A(:,1),2);
Q(:,1)=A(:,1)/R(1,1);
for j=2:n
    temp=Q;%找个地方储存
    Q(:,j)=A(:,j);
    R(:,j)=Q'*A(:,j);
    Q(:,j)=Q(:,j)-temp*R(:,j);%更新Q的第j列
    R(j,j)=norm(Q(:,j),2);
    Q(:,j)=Q(:,j)/R(j,j);  
end



%%% 用 MGS 方法实现 QR 分解
function [Q,R] = QR_MGS(A)
%
% 你的代码
%
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n);
R(1,1)=norm(A(:,1),2);
if R(1,1)~=0
    Q(:,1)=A(:,1)/R(1,1);
end

for j=2:n
    temp=Q;%找个地方储存
    Q(:,j)=A(:,j);
    R(:,j)=Q'*A(:,j);
    Q(:,j)=Q(:,j)-temp*R(:,j);%更新Q的第j列
    R(j,j)=norm(Q(:,j),2);
    if R(j,j)~=0
        Q(:,j)=Q(:,j)/R(j,j);  
    end
    
end


%%% 用 Householder 方法实现 QR 分解
function [Q,R] = QR_House(A)
%
% 你的代码
%
[m,n]=size(A);
Q=eye(m);
for k = 1:n
    x=A(k:m,k);
    [b,v]=House(x);
    A(k:m,k:n)= A(k:m,k:n)-b*v*(v'*A(k:m,k:n));
    Q(:,k:m)=Q(:,k:m)-b*(Q(:,k:m)*v)*v'  ;  
end
R=A;

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
