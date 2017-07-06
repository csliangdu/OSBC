function [weight,important_weight]=Chan_Attribute_weight(Data,beta)
% Author Cfy
% Creative Date 2008-11-24
% Modify Date 2009-05-07
% 参数说明 ：Data 为输入的数据集，beta是权重指标
% Function : 计算每一个属性在该数据集中的重要性，主要针对符号型数据
% 这个权重的计算利用Pattern Recognition 2004, (37): 943-952
[row,column]=size(Data);
%计算每一列对应的熵
entropy(1,column)=0;
important_weight(1,column)=0;
weight(1,column)=0;

% 统计每一列下属性值是否相同，判断有没有在某一列上属性值完全相同。
for i=1:column
    attributevalue=unique(Data(:,i));
    if length(attributevalue)==1
        % 在第i列上属性值相同，则为1        
        entropy(i)=1;
    else
        entropy(i)=0; 
    end;
end;
if length(find(entropy==0))==column
    %说明在每一列上属性值都是不同的，则计算每一列下的权重
    [Mode]=Find_Mode(Data); % 得到当前数据的mode
    for i=1:column
        % 得到当前列中Data中的所有对象到mode之间的距离
        i_column_value=row-length(find(Data(:,i)==Mode(:,i)));
        temp=0;
        for j=1:column
             j_column_value=row-length(find(Data(:,j)==Mode(:,j)));
             temp=temp+(i_column_value/j_column_value)^(1/(beta-1));
        end;
      weight(i)=1/temp;
    end;
else
    for i=1:column
        attributevalue=unique(Data(:,i));
        if length(attributevalue)==1
            weight(i)=1/length(find(entropy==1));
        else
             weight(i)=0;
        end;
    end;
end

% 如果某一列属性的权重大于均值，则认为该属性是重要的属性。
for i=1:column
     if weight(i)>=(1/column)
         important_weight(i)=1;    
     end;
 end;

    


 