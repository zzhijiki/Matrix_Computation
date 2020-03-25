R=[1,2,3,4,5;0,5,6,7,8;0,0,8,9,10;0,0,0,10,11;0,0,0,0,4];
v=[1,2,3,4,5]';
u=[4,5,6,7,8]';
u*v';
A=R+u*v'

n=length(v);
P=eye(n);
for i = n:-1:3
    G=eye(n);
    G1=givens(A(i-1,1),A(i,1));
    G(i-1:i,i-1:i)=G1;
    A=G*A;
    P=P*G';
end    
A   %验证是否为上HESS
P*A %验证回到了原来的矩阵


Q=eye(n);
for k=1:n
    for i=k+1:n
    G=eye(n);
    G1 = givens(A(k,k), A(i,k));

    B=[A(k,k:n);A(i,k:n)];
    B=G1*B;
    A(k,k:n)=B(1,:);
    A(i,k:n)=B(2,:);
    
    C=[Q(1:n,k),Q(1:n,i)];
    C=C*G1';
    Q(1:n,k)=C(:,1);
    Q(1:n,i)=C(:,2);
    
    end
end
A  %验证是否为上三角
P*Q*A %验证回到了原来的矩阵




