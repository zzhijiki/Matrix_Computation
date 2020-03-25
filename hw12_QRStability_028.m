%
% ���ѧ�ź�����
% 
% �����ַ���ʵ�� QR �ֽ⣬���������ȶ���
%
% �뽫�ļ����ͺ������е� xxx ��Ϊ���ѧ�������λ
%
function hw12_QRStability_028

m = 400; 
n = 300;

for cnd = 10.^[0:6]  % ������
    for k = 1 : 10   % ÿ������������ 10 ��
        A = randn(m,n);
        [U,S,V] = svd(A);  % ��������������Ͼ���
        
        % singular values range from 1 to cnd, with uniformly distributed logarithms
        sd = [1,cnd, exp(rand(1,n-2)*log(cnd))];
        A = U(:,1:n)*diag(sd)*V';   % ����ָ������ֵ���������
        
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
    
    %%% ������Խ��
    fprintf('\ncnd=%d\n',cnd);
    fprintf('   GS: relerr = %.4e, err_orth=%.4e\n', max(err_GS), max(err_GS_orth))
    fprintf('  MGS: relerr = %.4e, err_orth=%.4e\n', max(err_MGS),max(err_MGS_orth))
    fprintf('House: relerr = %.4e, err_orth=%.4e\n', max(err_H),  max(err_H_orth))
end

%%% �� GS ����ʵ�� QR �ֽ�
function [Q,R] = QR_GS(A)
%
% ��Ĵ���
%
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n);
R(1,1)=norm(A(:,1),2);
Q(:,1)=A(:,1)/R(1,1);
for j=2:n
    temp=Q;%�Ҹ��ط�����
    Q(:,j)=A(:,j);
    R(:,j)=Q'*A(:,j);
    Q(:,j)=Q(:,j)-temp*R(:,j);%����Q�ĵ�j��
    R(j,j)=norm(Q(:,j),2);
    Q(:,j)=Q(:,j)/R(j,j);  
end



%%% �� MGS ����ʵ�� QR �ֽ�
function [Q,R] = QR_MGS(A)
%
% ��Ĵ���
%
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(n);
R(1,1)=norm(A(:,1),2);
if R(1,1)~=0
    Q(:,1)=A(:,1)/R(1,1);
end

for j=2:n
    temp=Q;%�Ҹ��ط�����
    Q(:,j)=A(:,j);
    R(:,j)=Q'*A(:,j);
    Q(:,j)=Q(:,j)-temp*R(:,j);%����Q�ĵ�j��
    R(j,j)=norm(Q(:,j),2);
    if R(j,j)~=0
        Q(:,j)=Q(:,j)/R(j,j);  
    end
    
end


%%% �� Householder ����ʵ�� QR �ֽ�
function [Q,R] = QR_House(A)
%
% ��Ĵ���
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

function [beta,v]=House(x)  %������վ�ϸ���householder����
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
