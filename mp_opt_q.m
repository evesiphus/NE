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
mu = 300;

% number of sources
m = 3;

% number of players at each node
%N_m = sort(round(rand(1,m)*1000));
N = sort(randperm(3)*4,'descend');
n = sum(N);
% IP loss probability
Q=0.1:0.1:1;
p = 1-q;
% initinalization
TR_opt=zeros(1,length(Q));
U_opt=zeros(1,m);
V_opt=zeros(1,m);
U_all=[];
V_all=[];
TR_all=[];
for k = 1:length(Q)
    q=0.8;
   % q = (k-1).*0.1;
    p = 1-q;
    
    % initialization
    
    U=zeros(1,m);
    V=zeros(1,m);
    TR1=0;
    TR= 0;
    for ii=1:1:m
        V(1:ii)=0;
        i2=ii+1;
        U(i2:m)= N(i2:m);
        max_B =  sum(N(1:ii));
        
        for B=1:1:max_B
            sum_u_rest =max_B-B;
            x2=mod(sum_u_rest,2);
            x3=mod(sum_u_rest,3);
            x4=mod(sum_u_rest,4);
            % step (a)
            if ii==1
                U(ii)=sum_u_rest;
            elseif ii==2 && x2==0
                U(1)=sum_u_rest./2;
                U(2)=sum_u_rest./2;
            elseif ii==2 && x2==1
                U(1)=(sum_u_rest-1)./2+1;
                U(2)=(sum_u_rest-1)./2;
            elseif ii==3 && x3==0
                U(1)=(sum_u_rest)./3;
                U(2)=U(1);
                U(3)=U(1);
            elseif ii==3 && x3==1
                U(1)=(sum_u_rest-x3)./3+1;
                U(2)=(sum_u_rest-x3)./3;
                U(3)=(sum_u_rest-x3)./3;
            elseif ii==3 && x3==2
                U(1)=(sum_u_rest-x3)./3+1;
                U(2)=(sum_u_rest-x3)./3+1;
                U(3)=(sum_u_rest-x3)./3;
            elseif ii==4 && x4==1
                U(1)=(sum_u_rest-x4)./4+1;
                U(2)=(sum_u_rest-x4)./4;
                U(3)= U(2);
                U(4)= U(2);
            elseif ii==4 && x4==2
                U(1)=(sum_u_rest-x4)./4+1;
                U(2)=(sum_u_rest-x4)./4+1;
                U(3)= (sum_u_rest-x4)./4;
                U(4)= U(3);
            elseif ii==4 && x4==3
                U(1)=(sum_u_rest-x4)./4+1;
                U(2)=(sum_u_rest-x4)./4+1;
                U(3)=(sum_u_rest-x4)./4+1;
                U(4)=(sum_u_rest-x4)./4;
            end
            % step(b)
            jj=m-ii;
          if jj==1
              V(end)=B;
          elseif  jj==2
              if B<2
              V(end)=B;
              V(end-1)=0;
              else
               syms x y
              [x,y]=solve(N(end)+x*p==N(end-1)+y*p,x+y==B);
              V(end)=fix(x);
              V(end-1)=B-V(end);
              end
          elseif  jj==3 &&B==1
              V(end)=1;V(end-1)=0;V(end-2)=0;
          elseif  jj==3 &&B==2
              V(end)=1;V(end-1)=1;V(end-2)=0;    
          elseif  jj==3 && B>=3
              syms x y z
              [x,y,z]=solve(N(end)+x*p==N(end-1)+y*p,x+y+z==B,N(end)+x*p==N(end-2)+z*p);
              V(end)=fix(x);
              V(end-1)=fix(y);
              V(end-2)=B-V(end)-V(end-1);
          end
   end
              
            % step c
            if m==3
                TR =mu*(m-mu.*sum(1./(U(1).*phi+V(1).*phi.*p+mu)+ 1./(U(2).*phi+V(2).*phi.*p+mu)+1./(U(3).*phi+V(3).*phi.*p+mu)));
            elseif m==4
                TR =mu*(m-mu.*sum(1./(U(1).*phi+V(1).*phi.*p+mu)+1./(U(2).*phi+V(2).*phi.*p+mu)+1./(U(3).*phi+V(3).*phi.*p+mu)+1./(U(4).*phi+V(4).*phi.*p+mu)));
            end
            U_all=[U_all;U];
            V_all=[V_all;V];
            TR_all=[TR_all;TR];
            if TR>TR1
                TR1=TR;
                U_opt=U; V_opt=V;
            end
        end
   % end
    TR_opt(k)=TR1;
    TEXT=[U_all,V_all,TR_all]
end

figure(1)
plot(Q,TR_opt);
grid on;


