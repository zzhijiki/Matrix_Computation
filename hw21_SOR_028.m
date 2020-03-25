%
% SOR
%

function hw21_SOR_028

% 请将文件名和函数名中的 xxx 改为你的学号最后三位
%

N = 16; n = N*N;
a = 0; b = 1;
hx = (b-a)/(N+1);
hy = hx;
tol = 1e-6;
IterMax = 200;
u0 = zeros(N*N,1);   % initial guess

% coefficient matrix : A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%

% A=zeros(N*N)+diag((hx/2-1)*ones(N*N-1,1),1)+diag((hx/2-1)*ones(N*(N-1),1),N)+diag((4+hx*hx)*ones(N*N,1))+diag((-hx/2-1)*ones(N*N-1,1),-1)+diag((-hx/2-1)*ones(N*(N-1),1),-N);
% A(N+1:N:N*(N-1)+1,N:N:N*(N-1))=0;
% A(N:N:N*(N-1),N+1:N:N*(N-1)+1)=0

v1=ones(1,N);   %  Tn 的 diag
v2=ones(1,N-1);  %  Tn 的 次 diag
T=diag(2*v1)+diag(-1*v2,1)+diag(-1*v2,-1);   %  Tn
I=eye(N);  % I
D=diag(0.5*v2,1)+diag(-0.5*v2,-1);  %  D

A=kron(I,T)+kron(T,I)+hx.*kron(I,D)+hx.*kron(D,I)+(hx*hx).*kron(I,I);  %A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% right hand side


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  这里是你的代码  %%%%%

x=0:hx:1; x=x(2:N+1);
[X,Y]=meshgrid(x);
% f1 = (3+X.*Y).*(X+Y+X.*Y)-3*(X+Y).^2;  f1=f1(:);
% b=hx*hx*f1;
% norm_b=norm(b);
f=(3-2.*X).*(1-Y).*Y+(3-2.*Y).*(1-X).*X+(1-X).*X.*(1-Y).*Y;  f=f(:);
b=hx*hx*f;
norm_b=norm(b);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOR iteration
Omega = [1.5:0.05:1.9];
for k_omega = 1 : length(Omega)
    u = u0;
    omega = Omega(k_omega);
    for k = 1 : IterMax
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%  这里是你的代码  %%%%%
        
        for i=1:N*N
            ui_old=u(i);
            u(i)=1/A(i,i)*(b(i)-A(i,:)*u);
            u(i)=ui_old+omega*u(i);
        end

        relres = norm(b-A*u) / norm_b;
        if relres < tol
            break;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    Iter(k_omega) = k;
    fprintf('omega=%f, Iter=%d, relres=%.4e\n', omega, k, relres);
end
[itermin,idx] = min(Iter);
fprintf('最优参数为 %.2f, 对应的迭代步数为 %d\n',Omega(idx),itermin);

% plot the result
plot(Omega,Iter,'o-');
xlabel('\omega'); 
ylabel('iteration number');
legend([int2str(N),' x ', int2str(N)]);
title('SOR with different \omega');
% ===== End =====


