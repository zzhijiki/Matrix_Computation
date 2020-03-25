

function LU_tridiag
%ʹ���Լ���������
%A=[1,2,0;1,4,6;0,7,9]
%f=[2,3,4]'
%x=[3/4,5/8,-1/24]' ���ǵ�Ŀ��
a=[1,7];
b=[1,4,9];
c=[2,6];
f=[2,3,4];

n=length(b);
alpha=zeros(n,1);%��ʼ��һЩ����
beta=zeros(n-1,1);
y=zeros(n,1);
x=zeros(n,1);    %��ʼ��һЩ����


alpha(1)=b(1);
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