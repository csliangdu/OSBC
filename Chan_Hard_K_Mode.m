function cid=Chan_Hard_K_Mode(Data,K,beta)
% Author Cfy
% Date 2008-11-28
% Modify Date 2009-05-07
% ����˵�� ��Data Ϊ��������ݼ���InitialCenter�ǳ�ʼ��������,KΪҪ�����ĸ���,betaΪȨ��ָ��
% Function :ʵ��kmode�㷨��ģ���ӿռ������k-means�㷨��������빫ʽ���ü򵥵�0-1ƥ��
% importan_weight�Ǵ��ھ�ֵȨ������Ӧ������
% beta�Ƕ�Ӧ�Ĳ�������
InitialCenters=Centers(Data,K);
[n,d] = size(Data);
% ����cidΪ��������ʾ����
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
% ��ʼ����ÿһ�����Ե�Ȩ�س�ʼ��Ϊ1/d��
weight(K,d)=0;
important_weight(K,d)=0;
for j=1:K
    weight(j,:)=1/d;
    important_weight(j,:)=0;
end;

while 1== 1
%����ÿ�����ݵ��������ĵľ���
%Jobjection(iter+1)=0;
for i = 1:n
   %if bool==0 && ~isempty(find(InitialCenters==i))  cid(i)=find(InitialCenters==i); continue; end; %%�ж������������ͬ�Ķ������˳�ʼ�������ģ������ͻᵼ��ȱ��һ����
   temp=999999;
  for j=1:K
      dist= Subspace_Distance_of_Categorical(Data(i,:),ni(j,:),weight(j,:),beta);% ������������֮��ľ��뺯��   
      if dist<temp
          temp=dist;
          ind=j;
      end;   
  end
  cid(i) = ind;
end;
bool=1;
for i = 1:K
    %�ҵ�ÿһ����������ݣ��������ǵ�Mode����Ϊ�´μ���ľ������� K-mode�㷨
    ind = find(cid==i);
    %����ÿһ�������Զ�Ӧ��Ȩ��
    if ~isempty(ind)
        nj(i,:)=Find_Mode(Data(ind,:));
        [temp1,temp2]=Chan_Attribute_weight(Data(ind,:),beta);   
        weight(i,:)=temp1;
        important_weight(i,:)=temp2;
    end%nc(i,:) = mean(x(ind,:));
    % ͳ��ÿһ������ݸ���
    nr(i) = length(ind);
end
if ni==nj
    break 
else
    ni=nj;
 end
  iter = iter + 1;
end

