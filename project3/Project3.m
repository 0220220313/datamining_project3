clear all
for i=1:6
    data{i}=load(['graph_' int2str(i) '.txt']);
end
exdata=load('data3 .data');
data{7}=exdata(:,2:3)+1;
data{8}=data{7};
for i=1:8
    node=max(max(data{1,i}));
    T=spalloc(node,node,length(data{1,i}));
    for j=1:length(data{i})
        T(data{1,i}(j,1),data{1,i}(j,2))=1;
        if i==8
            T(data{1,i}(j,2),data{1,i}(j,1))=1;
        end
    end
    result{i}=[HITS(T,node,0.01,30),PageRank(T,node,0.0001,30)];
    if i<=5
        simrank{i}=SimRank(T,node,0.01,50,0.3);
    end
end
%% HITS
function result=HITS(T,node,delta,iter)
a_n=ones(node,1);
h_n=ones(node,1);
time=1;
while(time<iter)
    a_o=a_n;
    h_o=h_n;
    a_n=T'*h_o;
    h_n=T*a_n;
    a_n=a_n/norm(a_n);
    h_n=h_n/norm(h_n);
    if norm(a_n-a_o)+norm(h_n-h_o)<delta
        break;
    end
    time=time+1;
end
result=[a_n,SORT(a_n),h_n,SORT(h_n)];
end
%% PageRank
function result=PageRank(T,node,delta,iter)
D=0.15;
pr=ones(node,1)/node*D;
time=1;
while(time<iter)
    pr_o=pr;
    for i=1:node
        tmp=0;
        L=find(T(:,i));
        for j=1:length(L)
            tmp=tmp+pr_o(L(j))/length(find(T(L(j),:)));
        end
        pr(i)=D/node+(1-D)*tmp;
    end
    pr=pr/norm(pr);
    if norm(abs(pr-pr_o))<delta
        break;
    end
    time=time+1;
end
result=[pr/sum(pr)*100,SORT(pr)];
end
%% SimRank
function result=SimRank(T,node,delta,iter,C)
S=zeros(node,node);
S_o=S;
time=1;
for i=1:node
    L{i}=find(T(:,i));
end
while(time<iter)
    for i=1:node
        for j=i:node
            if(i==j)
                S(i,j)=1;
            elseif isempty(L{i})
                S(i,j)=0;
                break;
            elseif isempty(L{j})
                S(i,j)=0;
            else
                S(i,j)=C/length(L{i})/length(L{j});
                tmp=0;
                for i2=1:length(L{i})
                    for j2=1:length(L{j})
                        tmp=tmp+S(L{i}(i2),L{j}(j2));
                    end
                end
                S(i,j)=S(i,j)*tmp;
                S(j,i)=S(i,j);
            end
        end
    end
        if norm(S_o-S)<=delta
            break;
        end
    S_o=S;
    
    time=time+1;
end
result=[S,(1:node)'];
end
%% sort
function result=SORT(v)
[a,a_s]=sort(v);
result=zeros(length(v),1);
result(a_s(end))=1;
for i=length(v)-1:-1:1
    result(a_s(i))=length(v)+1-i;
    if a(i)==a(i+1)
        result(a_s(i))=result(a_s(i+1));
    end
end
end








