function [Mode]=Find_Mode(Data)
% Author Cfy
% Data 2008-02-13
% ����˵�� ��Data Ϊ��������ݼ�
% Function : �õ�һ�����ݼ��ϵ�mode
[row,column]=size(Data);
Mode=[];
for i=1:column
    ColumnValue=unique(Data(:,i));% �õ�ÿһ���в�ͬ������ֵ
    [Colrow, ColCol]=size(ColumnValue) ;% �õ���ǰ����ֵ��ά��
    ColumnValueSum=0; 
    dataset=[];
    for j=1:Colrow
       [Rowresult,Colresult]=size(find(Data(:,i)==ColumnValue(j,1)));
       %ColumnValueSum=Rowresult; % ��ʾ��j������ֵ�ڵ�i�г��ֵĴ���
       %ColumnValue(j,2)=ColumnValueSum  % �õ�ÿһ������ֵ���ֵĴ���
       dataset=cat(1,dataset,[ColumnValue(j,1),Rowresult]);
    end;
    [Maxrow,MaxValue]=max(dataset(:,2));
    Mode(1,i)=ColumnValue(MaxValue,1);
end;
    
       
        
    
    