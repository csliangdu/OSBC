function InitialCenters=Centers( data,k )

n=size(data,1);
cid=zeros(1,k);
for i=1:k
    %cid(i,:)=s(i,:); %直接取前三个元祖作为聚类中心
    mm=i*floor(n/k)-floor(rand(1,1)*(n/k));
    if mm==0
        mm=1;
    end
    cid(i)=mm;
end
InitialCenters=cid;
end

