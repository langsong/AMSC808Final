clear,clc;

%% Read data
filename = "lastfm_asia_edges.csv";
T = readtable(filename);
edgedata = table2array(T);

%% Problem 2
% s = edgedata(:,1)+1;
% t = edgedata(:,2)+1;
% G = graph(s,t);
% [bins, binsizes] = conncomp(G);

s = edgedata(:,1);
t = edgedata(:,2);
%G = graph(s,t);
n = max(t)+1;

adj = zeros(n); %initialize adjacency matrix
for j=1:length(s)
    adj(s(j)+1, t(j)+1)=1;
    adj(t(j)+1, s(j)+1)=1;
end

comp = DFS(adj);
disp(['Size of connected components are ', num2str(comp)])

%% Problem 3 (plot in problem 4 code)
deg = sum(adj,2); %degree for each vertex
degDist = tabulate(deg); %degree distribution
k = degDist(:,1);
pk = degDist(:,3)./100;

bin_k = zeros(10,1);
bin_pk = zeros(10,1);
product = k.*pk;

iteration = 0;
while true
    kmin = 2^(iteration);
    kmax = 2^(iteration+1);
    if kmax<= max(k)
        bin_k(iteration+1) = (sum(product(kmin:(kmax-1))))/(sum(pk(kmin:(kmax-1))));
        bin_pk(iteration+1) = (sum(pk(kmin:(kmax-1))))/2^iteration;
        iteration=iteration+1;
    else
%         bin_k(iteration+1) = (sum(product(kmin:end)))/(sum(pk(kmin:end)));
%         bin_pk(iteration+1) = (sum(pk(kmin:end)))/2^iteration;
        break
    end
end

%% Problem 4
lk = log(bin_k);
lpk = log(bin_pk); 
A = zeros(size(bin_k,1), 3);
A(:,1) = bin_k*(-1);
A(:,2) = lk*(-1);
A(:,3) = ones(size(bin_k,1),1);
A(any(isinf(A),2),:) = [];
lpk(any(isinf(lpk),2),:) = [];
b = A\lpk;
alpha = b(1);
tao = b(2);
C = exp(b(3));

sum_p = 0;
for i = 1:5000
    sum_p = sum_p + C*exp(-alpha*i)*(i^(-1*tao));
end
C = C/sum_p; %normalize C

%est_pk = C*exp(-alpha*bin_k).*(bin_k.^(-1*tao));
est_pk = C*exp(-alpha*k).*(k.^(-1*tao));

figure(1);
loglog(k,pk,'.', 'DisplayName','loglog points', 'MarkerSize',12)
hold on;
loglog(bin_k,bin_pk,'x','DisplayName','log-binning points','MarkerSize',12)
hold on;
%plot(bin_k,est_pk,'o','DisplayName','estimated degree distribution','MarkerSize',12)
plot(k,est_pk,'o','DisplayName','estimated degree distribution','MarkerSize',5)
xlabel('log(k)')
ylabel('log(pk)')
title('Loglog plot and log-binning plot of degree distribution')
legend()

%% Problem 5
s = s+1;
t = t+1;
G = graph(s,t);
d = distances(G);
ave_length = mean(d(~logical(eye(size(d)))));
disp(['The average shortest-path length in the network is ', num2str(ave_length)]);

z1 = 0; 
z2 = 0;
z3 = 0;
for i = 1:5000
    z1 = z1 + i*C*exp(-alpha*i)*(i^(-1*tao));
    z2 = z2 + (i^2)*C*exp(-alpha*i)*(i^(-1*tao));
    z3 = z3 + i*(i-1)*C*exp(-alpha*i)*(i^(-1*tao));
end
z2 = z2 - z1;
l = 1 + (log(n/z1))/(log(z2/z1));
disp(['The average shortest-path length for the random graph is ', num2str(l)]);

%% Problem 6
num_triangle = trace(adj^3)/6;

num_triple = 0;
for i = 1:length(adj)
    neighbors = find(adj(i,:)>0);
    if length(neighbors)<2; continue; end  %no triples
    num_triple = num_triple + nchoosek(length(neighbors),2);
end
num_triple = num_triple-2*num_triangle;
cluster_coef = 3*num_triangle/(num_triple+2*num_triangle);
cl_random = (1/n)*(z2^2)/(z1^3);
disp(['The coefficient for the actual network is ', num2str(cluster_coef)]);
disp(['The coefficient for the random network is ', num2str(cl_random)]);

%% Problem 8
product2 = product.*(k-1);
Tc = sum(product)/sum(product2);
disp(['Critical transmissibility for the actual network is ', num2str(Tc)]);

%% Problem 9
T = 0.4; 
perm = randperm(nnz(adj)/2);
num_edge = T*nnz(adj)/2;
count = 0;
%new adjacency matrix
for i = 1:length(adj)
    for j = i:length(adj)
        if(adj(i,j) ~= 0)
            count = count+1;
            adj(i,j) = (perm(count) <= num_edge);
        end
    end
end

for i = 1:length(adj)
    for j = 1:i
        adj(i,j) = adj(j,i);
    end
end

infected = zeros(10,1); %number of infections in 10 time steps
list = randperm(n);
start = list(1); 
I = zeros(n,1);
R = zeros(n,1);
I(start) = 1; 

for iter = 1:10
    cur = find(I); 
    %disp(cur);
    I(cur) = 0; 
    R(cur) = 1;
    s = length(cur);
    for j = 1:s
        neighbors = find(adj(cur(j),:)>0);
        nneigh = length(neighbors);
        for g = 1:nneigh
            if R(neighbors(g)) == 0
                I(neighbors(g)) = 1;
            end
        end
    end
    infected(iter) = sum(I);
end

timesteps = 1:10;
figure(2);
plot(timesteps,infected,'.','MarkerSize',15)
%plot(timesteps,infected)
xlabel('time step')
ylabel('number of infected')
title('Number of sick nodes versus time step')

figure(3);
%plot(timesteps,infected,'.','MarkerSize',15)
plot(timesteps,infected)
xlabel('time step')
ylabel('number of infected')
title('Number of sick nodes versus time step')

disp(['The fraction of nodes that are affected by the epidemic is ', num2str(sum(infected)/n)]);
%% Problem 10
Tc2 = z1/z3;
