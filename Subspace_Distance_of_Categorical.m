function Distance=Subspace_Distance_of_Categorical(Data1,Data2,weight,beta)
% Author Fuyuan Cao
% Date : 2008-11-24
% Modify Date 2009-05-07
% Function: ������������֮��ľ��뿼�����Ե���Ҫ�ԡ�
% beta�Ǳ仯�ģ�����Chan���㷨������Cao���㷨beta��1
[row,column]=size(Data1);
    Distance=0;
    for j=1:column
        if Data1(1,j)~=Data2(1,j)
            Distance=Distance+1*weight(j)^(beta);
        end;
    end;
        

