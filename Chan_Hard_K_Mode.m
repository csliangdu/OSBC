function cid=Chan_Hard_K_Mode(Data,K,beta)
% Author Cfy
% Date 2008-11-28
% Modify Date 2009-05-07
% 参数说明 ：Data 为输入的数据集，InitialCenter是初始聚类中心,K为要求聚类的个数,beta为权重指标
% Function :实现kmode算法，模仿子空间聚类中k-means算法，这个距离公式采用简单的0-1匹配
% importan_weight是大于均值权重所对应的属性
% beta是对应的参数因子
InitialCenters=Centers(Data,K);
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
% 初始化把每一个属性的权重初始化为1/d。
weight(K,d)=0;
important_weight(K,d)=0;
for j=1:K
    weight(j,:)=1/d;
    important_weight(j,:)=0;
end;

while 1== 1
%计算每个数据到聚类中心的距离
%Jobjection(iter+1)=0;
for i = 1:n
   %if bool==0 && ~isempty(find(InitialCenters==i))  cid(i)=find(InitialCenters==i); continue; end; %%判断如果有两个相同的对象都作了初始聚类中心，这样就会导致缺少一个类
   temp=999999;
  for j=1:K
      dist= Subspace_Distance_of_Categorical(Data(i,:),ni(j,:),weight(j,:),beta);% 调用两个对象之间的距离函数   
      if dist<temp
          temp=dist;
          ind=j;
      end;   
  end
  cid(i) = ind;
end;
bool=1;
for i = 1:K
    %找到每一类的所有数据，计算他们的Mode，作为下次计算的聚类中心 K-mode算法
    ind = find(cid==i);
    %计算每一类中属性对应的权重
    if ~isempty(ind)
        nj(i,:)=Find_Mode(Data(ind,:));
        [temp1,temp2]=Chan_Attribute_weight(Data(ind,:),beta);   
        weight(i,:)=temp1;
        important_weight(i,:)=temp2;
    end%nc(i,:) = mean(x(ind,:));
    % 统计每一类的数据个数
    nr(i) = length(ind);
end
if ni==nj
    break 
else
    ni=nj;
 end
  iter = iter + 1;
end

