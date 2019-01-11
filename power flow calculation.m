
pq=[1,2];%P,Q节点编号
pv=[3];
n=4;%节点数
m=4;%平衡节点编号

Y=[1.042093-8.242876i -0.588235+2.352941i 3.666667i -0.453858+1.891074i
    -0.588235+2.352941i 1.069005-4.727377i 0 -0.480769+2.403846i
    3.666667i 0 -3.333333i 0
    -0.453858+1.891074i -0.480769+2.403846i 0 0.934627-4.261590i];


%取电导和电纳
G=real(Y);B=imag(Y);
em=1.05;fm=0;%平衡节点电压

%初始化解
ve=[1 1 1.1];
vf=[0 0 0];
ps=[-0.3 -0.55 0.5 0];
qsvs=[-0.18 -0.13 1.21 0];%求解条件

f=ones(1,2*n);
maxi=1;
k=0;
while(maxi>0.1^5)
  k=k+1;
ue=[ve(1:m-1) em ve(m:n-1)];
uf=[vf(1:m-1) fm vf(m:n-1)];
%对PQ节点求f
for i=[m,pq]
    f(2*i-1)=ps(i)-ue(i)*(ue*G(:,i)-uf*B(:,i))-uf(i)*(uf*G(:,i)+ue*B(:,i));
    f(2*i)=qsvs(i)-uf(i)*(ue*G(:,i)-uf*B(:,i))+ue(i)*(uf*G(:,i)+ue*B(:,i));
end
%对pv节点求f
for i=pv
    f(2*i-1)=ps(i)-ue(i)*(ue*G(:,i)-uf*B(:,i))-uf(i)*(uf*G(:,i)+ue*B(:,i));
    f(2*i)=qsvs(i)-uf(i)^2-ue(i)^2;
end
ff=f;
ff(2*m-1:2*m)=[];
maxi=max(abs(ff));
    
%求雅可比矩阵   
J=zeros(2*n,2*n);
for i=pq
    for j=[1:m-1,m+1:n]
    if(i==j)
        J(2*i-1,2*j-1)=-(ue*G(:,i)-uf*B(:,i))-G(i,i)*ue(i)-B(i,i)*uf(i);%dP/de
        J(2*i-1,2*j)=-(uf*G(:,i)+ue*B(:,i))+B(i,i)*ue(i)-G(i,i)*uf(i);%dP/df
        J(2*i,2*j-1)=(uf*G(:,i)+ue*B(:,i))+B(i,i)*ue(i)-G(i,i)*uf(i);%dQ/de
        J(2*i,2*j)=-(ue*G(:,i)-uf*B(:,i))+G(i,i)*ue(i)+B(i,i)*uf(i);%dQ/df
    else
        J(2*i-1,2*j-1)=-G(i,j)*ue(i)-B(i,j)*uf(i);
        J(2*i,2*j)=-J(2*i-1,2*j-1);%dP/de=-dQ/df
        J(2*i-1,2*j)=B(i,j)*ue(i)-G(i,j)*uf(i);
        J(2*i,2*j-1)=J(2*i-1,2*j);%dP/df=dQ/de
    end
    end
end
for i=pv
    for j=[1:m-1,m+1:n]
    if(i==j)
        J(2*i-1,2*j-1)=-(ue*G(:,i)-uf*B(:,i))-G(i,i)*ue(i)-B(i,i)*uf(i);%dP/de
        J(2*i-1,2*j)=-(uf*G(:,i)+ue*B(:,i))+B(i,i)*ue(i)-G(i,i)*uf(i);%dP/df
        J(2*i,2*j-1)=-2*ue(i);%dV/de
        J(2*i,2*j)=-2*uf(i);%dV/df
    else
        J(2*i-1,2*j-1)=0;J(2*i,2*j)=0;
        J(2*i-1,2*j)=0;J(2*i,2*j-1)=0;
    end
    end
end
JJ=J;
JJ(2*m-1:2*m,:)=[];
JJ(:,2*m-1:2*m)=[];

detx=-JJ\ff';

%更新电压纵分量和横分量
for i=1:n-1
ve(i)=ve(i)+detx(2*i-1);
vf(i)=vf(i)+detx(2*i);
end
end