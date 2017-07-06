function new_class=k_modes_ndm(data,class)
t0=clock;
[r,cc]=size(data);
center=Centers( data,class );
k=length(center);
modes=zeros(k,cc);
FSS=[];
NEWFSS=0;
OLDFSS=0;
clen=zeros(1,cc);
for i=1:cc
    column=unique(data(:,i));
    clen(i)=length(column);
    for j=1:clen(i)
        attr(j,i)=column(j);
    end
end
[ar,ac]=size(attr);

for i=1:k
    modes(i,:)=data(center(i),:);
end
old_mode=zeros(k,ar,cc);
for i=1:k
    for j=1:cc
        w=find(attr(:,j)==modes(i,j));
        x=w(1);
        old_mode(i,x,j)=1;
    end
end
    
t=0;
new_class=zeros(1,r);
old_class=zeros(1,r);
while 1==1
    if t>40
        break;
    end
    t=t+1;
    NEWFSS=0;
    new_mode=zeros(k,ar,cc);
    classsum=zeros(1,k);
    for i=1:r
        distance=zeros(1,k);
        for j=1:k
            dis=0;
            for h=1:cc
                for l=1:clen(h)
                    if attr(l,h)==data(i,h)
                        dis=dis+(1-old_mode(j,l,h))^2;
                    else
                        dis=dis+old_mode(j,l,h)^2;
                    end
                end
            end
            distance(j)=distance(j)+dis;
        end
        [MinValue,MinRow]=min(distance);
        old_class(i)=MinRow;
        NEWFSS=NEWFSS+MinValue;
        classsum(MinRow)=classsum(MinRow)+1;
        for e=1:cc
            w=find(attr(:,e)==data(i,e));
            x=w(1);
            new_mode(MinRow,x,e)=new_mode(MinRow,x,e)+1;
        end
    end
   
    for i=1:k
        if classsum(i)~=0
            old_mode(i,:,:)=new_mode(i,:,:)/classsum(i);
        
        end
    end
    if OLDFSS==NEWFSS 
        FSS(t)=NEWFSS;
        break;
    else
        OLDFSS=NEWFSS;
          new_class=old_class;
          FSS(t)=NEWFSS;
          %oldnum=classsum;
    end
%     if isequal(old_class,new_class)
%         break;
%     else
%         new_class=old_class;
%     end
end
times=etime(clock,t0);
end