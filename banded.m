n=6; %���ɴ�״����ķ����������еĶ����� �����Խ��� �� �ζԽ���
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