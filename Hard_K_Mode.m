function cid=Hard_K_Mode(Data,K)
InitialCenters=Centers(Data,K);
%InitialCenters=inc(xu,:);
% Author Cfy
% Data 2008-09-19
% 参数说明 ：Data 为输入的数据集，InitialCenter是初始聚类中心,K为要求聚类的个数
% Function :实现kmode算法，模仿K-means算法
[n,d] = size(Data);
% 设置cid为分类结果显示矩阵
cid = zeros(1,n); 
% Make this different to get the loop started.
oldcid = ones(1,n);
% The number in each cluster.
nr = zeros(1,K); 
% Set up maximum number of iterations.
maxgn= 100;
iter = 0;
InitialCenter=Data(InitialCenters,:);
ni=InitialCenter;
bool=0;
while 1 == 1
%计算每个数据到聚类中心的距离
%Jobjection(iter+1)=0;
for i = 1:n
   if bool==0 && ~isempty(find(InitialCenters==i))  cid(i)=find(InitialCenters==i); continue; end; %%判断如果有两个相同的对象都作了初始聚类中心，这样就会导致缺少一个类
   dist= Distance_of_Categorical(repmat(Data(i,:),K,1),ni);% 调用两个对象之间的距离函数
   [m,ind] = min(dist(:,2)); % 将当前聚类结果存入cid中
    %Jobjection(iter+1)=Jobjection(iter+1)+m;
   cid(i) = ind;
end
bool=1;
nj=ni;
for i = 1:K
%找到每一类的所有数据，计算他们的Mode，作为下次计算的聚类中心 K-mode算法
ind = find(cid==i);
if ~isempty(ind)
nj(i,:)=Find_Mode(Data(ind,:));
end%nc(i,:) = mean(x(ind,:));
% 统计每一类的数据个数
nr(i) = length(ind);
end
 if ni==nj break 
 else
       ni=nj;
 end
  iter = iter + 1;
end

