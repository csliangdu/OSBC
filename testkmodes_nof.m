function testkmodes_nof(s,k,ax)
tic;
diedai=100;
for xu=1:diedai
[n,p]=size(s);
s1=s(:,1:p-1);
Ak=s(:,p);
cluster=K_modes_nof(s1,k,ax);
cluster=cluster';
Re=[cluster,Ak];
%--------------------求ARI和ac
[n,p]=size(Re);
A1=Re(:,1)';
A2=Re(:,2)';
cp=zeros(k,k);
for i=1:k
    p1=find(A1==i);
    for j=1:k
      c1=find(A2==j);
      cp(i,j)=length(intersect(p1,c1));
    end
end
cp(k+1,:)=sum(cp);
cp(:,k+1)=sum(cp,2);
%for j=1:k
%    [a,b]=size(find(A1==j));
%    p(j)=a;
%end
%for j=1:k
%    [a,b]=size(find(A2==j));
%    c(j)=a;
%end
%-----------------------------------------------------------
for i=1:k  %计算ac
    a(i)=max(cp(i,1:k));
end
ac(xu)=sum(a)/n
%---------------------------------------------------------------
%-----------------------------------------------
%计算ARI
cp1=cp(1:k,1:k);
r0=sum(sum((cp1.*(cp1-1))./2));
cp2=cp(1:k,1+k)';
r1=sum((cp2.*(cp2-1))./2);
cp3=cp(1+k,1:k);
r2=sum((cp3.*(cp3-1))./2);
r3=(2*r1*r2)/(n*(n-1));
ARI(xu)=(r0-r3)/(0.5*(r1+r2)-r3)
%--------------------------------------------
%---------------------------计算纯度,精度
pre(xu)=(sum(max(cp(1:k,1:k),[],2)./cp(1:k,k+1)))/k
re(xu)=(sum(max(cp(1:k,1:k))./cp(k+1,1:k)))/k
%-----------------------------------
end
AC(1)=sum(ac)/diedai;
AC(2)=max(ac);
AC(3)=min(ac);
AC(4)=sum((ac-AC(1)).^2)/diedai
ari(1)=sum(ARI)/diedai;
ari(2)=max(ARI);
ari(3)=min(ARI);
ari(4)=sum((ARI-ari(1)).^2)/diedai
Pre(1)=sum(pre)/diedai;
Pre(2)=max(pre);
Pre(3)=min(pre);
Pre(4)=sum((Pre(1)-pre).^2)/diedai
RE(1)=sum(re)/diedai;
RE(2)=max(re);
RE(3)=min(re);
RE(4)=sum((RE(1)-re).^2)/diedai
%plot(convergence);figure(gcf);
toc

end

