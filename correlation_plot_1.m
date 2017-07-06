function correlation_plot_1(data)
[similariyt,similaritycos,similarityEu]=matrix(data);
quxian1(similariyt,similaritycos,similarityEu);
end

function [similariyt,similaritycos,similarityEu]=matrix(A)
[n,d]=size(A);
similarityper=zeros(n);
for i=1:n 
    xi=A(i,:);
    similarityper(i,i)=0;
    for j=i+1:n
        xj=A(j,:);
       similarity=(xi~=xj);
       similarityper(i,j)=sum(similarity)/(d-1);  
    end
end
v=diag(similarityper);
similariyt=similarityper+diag(v)*(-1)+similarityper';
maxx=max(max(similariyt));
minn=min(min(similariyt));
similariyt=similariyt./(maxx-minn);
[similaritycos]=matrix_cos(similariyt);
[similarityEu]=matrix_Eu(similariyt);
end


function [similarityEu]=matrix_Eu(A)
n=size(A,1);
similarityper=zeros(n);
for i=1:n  %构造相似比矩阵
    xi=A(i,:);
    for j=i+1:n
        xj=A(j,:);
       similarityper(i,j)=(sqrt(dot((xi-xj),(xi-xj))));  
    end
end
v=diag(similarityper);
similarityEu=similarityper+diag(v)*(-1)+similarityper';
maxx=max(max(similarityEu));
minn=min(min(similarityEu));
similarityEu=similarityEu./(maxx-minn);
end


function [similaritycos]=matrix_cos(A)
n=size(A,1);
similarityper=zeros(n);
for i=1:n  %构造相似比矩阵
    xi=A(i,:);
    for j=i+1:n
        xj=A(j,:);
       similarityper(i,j)=2*sqrt(1-(dot(xi,xj)./(sqrt(dot(xi,xi)).*sqrt(dot(xj,xj)))));
    end
end
v=diag(similarityper);
similaritycos=similarityper+diag(v)*(-1)+similarityper';
maxx=max(max(similaritycos));
minn=min(min(similaritycos));
similaritycos=similaritycos./(maxx-minn);
end


function quxian1(similariyt,similaritycos,similarityEu)
vec=xiasanjiao(similariyt);
veccos=xiasanjiao(similaritycos);
veceu=xiasanjiao(similarityEu);
[sortvec,locat]=sort(vec);
sortveccos=veccos(locat);
sortveceu=veceu(locat);
xn=length(sortvec);
x=1:1:xn;
figure();
subplot(3,1,1);
plot(x,sortvec,'r');
title('0-1 distance on the original data');
subplot(3,1,2);
plot(x,sortveccos);
title('Cosine distance on the space structure');
subplot(3,1,3);
plot(x,sortveceu,'Color',[0 0.5 0]);
title('Euclidean distance on the space structure');
end


function [vector]=xiasanjiao(sim)
n=size(sim,1);
vector=[];
for i=1:n
    vector=[vector sim(i,i+1:end)];
end
end
