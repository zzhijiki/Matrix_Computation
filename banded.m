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