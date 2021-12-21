function componentsize=DFS(A)
% DFS Algorithm Modified from lecture notes
% Added extra steps and the output is a vector of individual component
% sizes
    n=size(A,1);
    componentsize=zeros(n,1);
    curr_size=0;
    color=zeros(n,1);
    parent=zeros(n,1);
    d=zeros(n,1);
    f=zeros(n,1);
    time=0;
    for i=1:n
        if color(i)==0
            [time,color,parent,d,f]=DFSVISIT(A,i,time,color,parent,d,f);
            curr_size=curr_size+1;
            if curr_size==1
                componentsize(curr_size)=sum(color);
            else
                componentsize(curr_size)=sum(color)-sum(componentsize(1:( ...
                    curr_size-1)));
            end
        end
    end
    componentsize=componentsize(1:curr_size); 
end

function [time,color,parent,d,f]=DFSVISIT(A,u,time,color,parent,d,f)
    time=time+1;
    d(u)=time;
    color(u)=-1;
    i=find(A(u,:)>0);
    for j=1:length(i)
        if color(i(j))==0
            parent(i(j))=u;
            [time,color,parent,d,f]=DFSVISIT(A,i(j),time,color,parent,d,f);
        end
    end
    color(u)=1;
    time=time+1;
    f(u)=time;
end