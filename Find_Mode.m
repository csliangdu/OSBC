function [Mode]=Find_Mode(Data)
% Author Cfy
% Data 2008-02-13
% 参数说明 ：Data 为输入的数据集
% Function : 得到一个数据集合的mode
[row,column]=size(Data);
Mode=[];
for i=1:column
    ColumnValue=unique(Data(:,i));% 得到每一列中不同的属性值
    [Colrow, ColCol]=size(ColumnValue) ;% 得到当前属性值得维数
    ColumnValueSum=0; 
    dataset=[];
    for j=1:Colrow
       [Rowresult,Colresult]=size(find(Data(:,i)==ColumnValue(j,1)));
       %ColumnValueSum=Rowresult; % 表示第j个属性值在第i列出现的次数
       %ColumnValue(j,2)=ColumnValueSum  % 得到每一个属性值出现的次数
       dataset=cat(1,dataset,[ColumnValue(j,1),Rowresult]);
    end;
    [Maxrow,MaxValue]=max(dataset(:,2));
    Mode(1,i)=ColumnValue(MaxValue,1);
end;
    
       
        
    
    