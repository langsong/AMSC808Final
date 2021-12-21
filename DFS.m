function comp = DFS(adj)
    n = size(adj,1);
    comp = zeros(n,1);
    cur_size = 0;
    color = zeros(n,1);
    parent = zeros(n,1);
    d = zeros(n,1);
    f = zeros(n,1);
    time = 0;
    for i = 1:n
        if color(i) == 0
            [time,color,parent,d,f] = DFSVISIT(adj,i,time,color,parent,d,f);
            cur_size = cur_size+1;
            if cur_size == 1
                comp(cur_size) = sum(color);
            else
                comp(cur_size) = sum(color)-sum(comp(1:(cur_size-1)));
            end
        end
    end
    comp = comp(1:cur_size); 
end

function [time,color,parent,d,f] = DFSVISIT(A,u,time,color,parent,d,f)
    time = time+1;
    d(u) = time;
    color(u) = -1;
    i = find(A(u,:)>0);
    for j = 1:length(i)
        if color(i(j)) == 0
            parent(i(j)) = u;
            [time,color,parent,d,f] = DFSVISIT(A,i(j),time,color,parent,d,f);
        end
    end
    color(u) = 1;
    time = time+1;
    f(u) = time;
end