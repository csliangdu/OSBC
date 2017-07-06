function correlation_plot_2(data)
[similariyt,similaritycos,similarityEu]=matrix(data);
[vec,veccos,veceu]=Untitled7(similariyt,similaritycos,similarityEu);
coplot(vec,veccos,veceu);
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
for i=1:n 
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
for i=1:n  
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

function [vec,veccos,veceu]=Untitled7(similariyt,similaritycos,similarityEu)
vec=xiasanjiao(similariyt);
veccos=xiasanjiao(similaritycos);
veceu=xiasanjiao(similarityEu);
end


function [vector]=xiasanjiao(sim)
n=size(sim,1);
vector=[];
for i=1:n
    vector=[vector sim(i,i+1:end)];
end
end

function coplot(vec,veccos,veceu)
figure();
subplot(1,2,1);
plot(vec,veccos,'.');
xlabel('0-1 distance on the original data');
ylabel('Cosine distance on the space structure');
subplot(1,2,2);
plot(vec,veceu,'.','Color',[0 0.5 0]);
xlabel('0-1 distance on the original data');
ylabel('Euclidean distance on the space structure');
end