function [weight,important_weight]=Chan_Attribute_weight(Data,beta)
% Author Cfy
% Creative Date 2008-11-24
% Modify Date 2009-05-07
% ����˵�� ��Data Ϊ��������ݼ���beta��Ȩ��ָ��
% Function : ����ÿһ�������ڸ����ݼ��е���Ҫ�ԣ���Ҫ��Է���������
% ���Ȩ�صļ�������Pattern Recognition 2004, (37): 943-952
[row,column]=size(Data);
%����ÿһ�ж�Ӧ����
entropy(1,column)=0;
important_weight(1,column)=0;
weight(1,column)=0;

% ͳ��ÿһ��������ֵ�Ƿ���ͬ���ж���û����ĳһ��������ֵ��ȫ��ͬ��
for i=1:column
    attributevalue=unique(Data(:,i));
    if length(attributevalue)==1
        % �ڵ�i��������ֵ��ͬ����Ϊ1        
        entropy(i)=1;
    else
        entropy(i)=0; 
    end;
end;
if length(find(entropy==0))==column
    %˵����ÿһ��������ֵ���ǲ�ͬ�ģ������ÿһ���µ�Ȩ��
    [Mode]=Find_Mode(Data); % �õ���ǰ���ݵ�mode
    for i=1:column
        % �õ���ǰ����Data�е����ж���mode֮��ľ���
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

% ���ĳһ�����Ե�Ȩ�ش��ھ�ֵ������Ϊ����������Ҫ�����ԡ�
for i=1:column
     if weight(i)>=(1/column)
         important_weight(i)=1;    
     end;
 end;

    


 