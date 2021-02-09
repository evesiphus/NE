%% wireless network games: optimal solution of multi sources with pure strategies (ALG)
% * <index.html *INDEX*>



%%
%clc;
%clear;
%close all;

%% parameter seeting

%traffic generation rate
phi = 1;

%traffic average service time
mu = 1;

% number of sources
m = 3;

% number of players at each node
%N_m = sort(round(rand(1,m)*1000));
%N = [50,5,5,5,5,5];
N=[12 8 4];
n = sum(N);

%traffic average service time
Mu = 0:1:10;

% initinalization
TR_opt=zeros(1,length(Mu));
U_opt=zeros(1,m);
V_opt=zeros(1,m);
U_all=[];
V_all=[];
TR_all=[];
for k = 1:length(Mu) %length(Q)=11, q=0,0.1,...,1
    mu=1.*(k-1);
    q = 0.5;
    p = 1-q;
    
    % initialization
    U=zeros(1,m);
    V=zeros(1,m);
    
    %initilize TR1 as the traffic when all players choose DP
    LL=zeros(1,m);
    for h=1:m
        LL(h)=mu./(N(h)*phi+mu);
    end
    TR1=mu*(m-sum(LL));
    U_opt=N; V_opt=zeros(1,m);
    
    for ii=1:1:m-1 % changed to m-1 !!!!!
        V(1:ii)=0;  %V(1)=...=V(ii)=0
        i2=ii+1;
        U(i2:m)= N(i2:m);
        max_B =  sum(N(1:ii));
        
        for B=1:1:max_B
            sum_u_rest =max_B-B;
            
            % step (a)
            U(1:ii)=N(1:ii);
            for t=1:1:B
                M=U(1:ii);
                index_set=find(M==max(M));
                ind=index_set(1);
                U(ind)=U(ind)-1;
            end
            
            %step (b)
            V(i2:m)=0;
            for t=1:1:B
                W=zeros(1,m);
                for j=i2:m
                    W(j)=N(j)*phi+V(j)*p*phi;
                end
                M=W(i2:m);
                index_set=find(M==min(M));
                ind=index_set(1);
                V(ii+ind)=V(ii+ind)+1;
            end
            
            
            % step c
            L=zeros(1,m);
            for h=1:m
                L(h)=mu./(U(h)*phi+V(h)*p*phi+mu);
            end
            TR=mu*(m-sum(L));
            
            if TR>TR1           %update
                TR1=TR;
                U_opt=U; V_opt=V;
            end
        end
    end
    
    TR_opt(k)=TR1;
    TEXT=[k,U_opt,V_opt,TR_opt(k)];
end

figure(10)
plot(Mu,TR_opt);
grid on;


