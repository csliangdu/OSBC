function SBC_1(data,k,t)
%--------------------------------------------------------------------------
%  input:
%       data:data contain label information, but label information only
%            used in evaluate process
%       k: classes parameter
%       t: running times
%
%  output:
%       t times cluster results was saved in cl.mat
%       t times evaluation index was saved in index.mat
%--------------------------------------------------------------------------
tstart=tic;
[n1,m1]=size(data);
cl=zeros(n1,t);
Ak=data(:,m1);
A=data(1:n1,1:m1-1);
[n,m]=size(A);
similarityper=zeros(n);
convergence=zeros(t,10);
tstarts=tic;
%--------------------------------------------------------------------------
%   create space structure 
for i=1:n  
    xi=A(i,:);
    similarityper(i,i)=1;
    for j=i+1:n
        xj=A(j,:);
       similarity=(xi==xj);
       similarityper(i,j)=sum(similarity)/m;  
    end
end
v=diag(similarityper);
s=similarityper+diag(v)*(-1)+similarityper';
%--------------------------------------------------------------------------
initialcenter=zeros(t,k);
for j=1:t
    for i=1:k
          mm=i*floor(n/k)-floor(rand(1,1)*(n/k));
          if mm==0
               mm=1;
          end
          initialcenter(j,i)=mm;
    end
end
mtime=toc(tstarts);
ac=zeros(1,t);
ARI=zeros(1,t);
pre=zeros(1,t);
re=zeros(1,t);
for xu=1:t
    fprintf('NUMBER OF CLUSTER TIMES: %i \n', xu);
[n,p]=size(s);
s(:,p+1)=0;
cid=s(initialcenter(xu,:),:);
flags=1;
times=1;
cishu=1;
while flags
    flags=0;
    times=times+1; 
    dist=zeros(n,k);
     for i=1:n
         for j=1:k
         dist(i,j)=sqrt(2*(1-(dot(s(i,:),cid(j,:))./(sqrt(dot(s(i,:),s(i,:))).*sqrt(dot(cid(j,:),cid(j,:)))))));
         end
         [~,y]=find(dist(i,:)==min(dist(i,:)));
         c=size(find(y(1)==s(i,p+1)),1);
         if c==0 
         flags=flags+1;    
         s(i,p+1)=y(1);
         else 
             continue;
         end
     end
     Csum=zeros(1,k);
     flag2=0;
     for j=1:k
         Asum=0;
         r=find(s(:,p+1)==j);
         if isempty(r)
             [~,lot]=find(dist(:,j)==max(dist(:,j)));
             cid(j,:)=s(lot(1),:);
             cid(j,end)=j;
             flag2=1;
         else
         cid(j,:)=mean(s(r,:),1);
         end
         for m=1:length(r)
             Asum=Asum+sqrt(2*(1-(dot(s(r(m),:),cid(j,:))./(sqrt(dot(s(r(m),:),s(r(m),:))).*sqrt(dot(cid(j,:),cid(j,:)))))));
         end
         Csum(1,j)=Asum;
     end
     if flag2==1
         continue;
     end
     sum(Csum(1,:));
     Csum2=sum(Csum(1,:));
     convergence(xu,cishu)=Csum2;
     cishu=cishu+1;
end
convergence(xu,cishu:10)=Csum2;
cluster=s(:,p+1);
cl(:,xu)=cluster;
Re=[cluster,Ak];

%--------------------compute AC,ARI,PRE,RE-------------------------
n=size(Re,1);
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
%------------------------------------------------------------------
a=zeros(1,k);
for i=1:k  
    a(i)=max(cp(i,1:k));
end
ac(xu)=sum(a)/n;
%-------------------------------------------------------------------

%------------------------------------------------------------------
cp1=cp(1:k,1:k);
r0=sum(sum((cp1.*(cp1-1))./2));
cp2=cp(1:k,1+k)';
r1=sum((cp2.*(cp2-1))./2);
cp3=cp(1+k,1:k);
r2=sum((cp3.*(cp3-1))./2);
r3=(2*r1*r2)/(n*(n-1));
ARI(xu)=(r0-r3)/(0.5*(r1+r2)-r3);
%------------------------------------------------------------------

%----------------------------------------------------------------
pre(xu)=(sum(max(cp(1:k,1:k),[],2)./cp(1:k,k+1)))/k;
re(xu)=(sum(max(cp(1:k,1:k))./cp(k+1,1:k)))/k;
%----------------------------------------------------------------
end
AC=zeros(1,4);
AC(1)=sum(ac)/t;
AC(2)=max(ac);
AC(3)=min(ac);
AC(4)=sum((ac-AC(1)).^2)/t;
ari=zeros(1,4);
ari(1)=sum(ARI)/t;
ari(2)=max(ARI);
ari(3)=min(ARI);
ari(4)=sum((ARI-ari(1)).^2)/t;
Pre=zeros(1,4);
Pre(1)=sum(pre)/t;
Pre(2)=max(pre);
Pre(3)=min(pre);
Pre(4)=sum((Pre(1)-pre).^2)/t;
RE=zeros(1,4);
RE(1)=sum(re)/t;
RE(2)=max(re);
RE(3)=min(re);
RE(4)=sum((RE(1)-re).^2)/t;
convergence=real(convergence)';
convergence=convergence(1:10,:);
plot(convergence);
set(gca,'ytick',[]);
xlabel('Iteration times');
save cl cl;

fprintf('\r Save cluster results in Current Folder named cl.mat.\r\n');
fprintf('index of AC:mean:  %4.4f    max: %4.4f    min: %4.4f    variance: %4.4f \r\n', AC(1),AC(2),AC(3),AC(4));
fprintf('index of ARI:mean: %4.4f    max: %4.4f    min: %4.4f    variance: %4.4f \r\n', ari(1),ari(2),ari(3),ari(4));
fprintf('index of PRE:mean: %4.4f    max: %4.4f    min: %4.4f    variance: %4.4f \r\n', Pre(1),Pre(2),Pre(3),Pre(4));
fprintf('index of RE:mean:  %4.4f    max: %4.4f    min: %4.4f    variance: %4.4f \r\n', RE(1),RE(2),RE(3),RE(4));
index=[AC;ari;Pre;RE];
save index index

fprintf('Save index result in Current Folder named index.mat.\r\n');

%-------------compute Average computational time-------------------------
alltime=toc(tstart);
averagetime=(alltime-mtime)/t+mtime;
fprintf('Average cost time is %f seconds, built matrix cost %f seconds.\n',averagetime,mtime);

end
