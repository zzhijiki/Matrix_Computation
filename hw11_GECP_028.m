%
% ���ѧ�ź�����
% 
% ȫ��Ԫ LU �ֽ�: P*A*Q = L*U
%
% �뽫�ļ����ͺ������е� xxx ��Ϊ���ѧ�������λ
%
function hw11_GECP_028  

n = 100;

A = randn(n); % ���� n ���������

%%% ȫ��Ԫ LU �ֽ�
[L,U,p,q] = GECP(A); 

%%% ���ԣ�����������  ||P*A*Q - L*U|| / ||A||�������� F-����
fprintf('relerr = %.4e\n', norm( A(p,q) - L*U, 'fro') / norm(A,'fro') )


%%% ȫ��Ԫ LU �ֽ�: P*A*Q = L*U
function [L,U,p,q] = GECP(A)
%
% ��Ĵ���
%
[n,n]=size(A);
p=1:n;  q=1:n;
for i = 1: n-1
    
    B=A(i:n,i:n);
    [a,I]=max(abs(B(:)));
    
    [k_row, k_col] = ind2sub(size(B),I);
    if a== 0
        erorr("Error:��%d��������ԪΪ0��\n", i)
    end
    k_row=k_row+i-1;  k_col=k_col+i-1;
    if k_row~=1
        temp=A(i,:);A(i,:)=A(k_row,:);A(k_row,:)=temp;
        temp=p(i);p(i)=p(k_row);p(k_row)=temp;
    end
    if k_col~=1
        temp=A(:,i);A(:,i)=A(:,k_col);A(:,k_col)=temp;
        temp=q(i);q(i)=q(k_col);q(k_col)=temp;
    end  

    A(i+1:n,i)=A(i+1:n,i)/A(i,i);
    A(i+1:n,i+1:n)=A(i+1:n,i+1:n)-A(i+1:n,i)*A(i,i+1:n);
end
U=triu(A);
L=eye(n)+tril(A,-1);


    
