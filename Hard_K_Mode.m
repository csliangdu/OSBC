function cid=Hard_K_Mode(Data,K)
InitialCenters=Centers(Data,K);
%InitialCenters=inc(xu,:);
% Author Cfy
% Data 2008-09-19
% ����˵�� ��Data Ϊ��������ݼ���InitialCenter�ǳ�ʼ��������,KΪҪ�����ĸ���
% Function :ʵ��kmode�㷨��ģ��K-means�㷨
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
while 1 == 1
%����ÿ�����ݵ��������ĵľ���
%Jobjection(iter+1)=0;
for i = 1:n
   if bool==0 && ~isempty(find(InitialCenters==i))  cid(i)=find(InitialCenters==i); continue; end; %%�ж������������ͬ�Ķ������˳�ʼ�������ģ������ͻᵼ��ȱ��һ����
   dist= Distance_of_Categorical(repmat(Data(i,:),K,1),ni);% ������������֮��ľ��뺯��
   [m,ind] = min(dist(:,2)); % ����ǰ����������cid��
    %Jobjection(iter+1)=Jobjection(iter+1)+m;
   cid(i) = ind;
end
bool=1;
nj=ni;
for i = 1:K
%�ҵ�ÿһ����������ݣ��������ǵ�Mode����Ϊ�´μ���ľ������� K-mode�㷨
ind = find(cid==i);
if ~isempty(ind)
nj(i,:)=Find_Mode(Data(ind,:));
end%nc(i,:) = mean(x(ind,:));
% ͳ��ÿһ������ݸ���
nr(i) = length(ind);
end
 if ni==nj break 
 else
       ni=nj;
 end
  iter = iter + 1;
end

