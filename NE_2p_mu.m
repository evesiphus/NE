%% wireless network games: Nash equilibrium of two sources with pure strategies
% * <index.html *INDEX*>
 
 
 
%%
%clc;
%clear;
%close all;


%% parameter seeting

%traffic generation rate
phi = 1;

%traffic average service time
Mu = 100:100:6000;

% number of players at each node
n_1 = 1000;
n_2 =100;

% IP loss probability
%Q = 0.1:0.1:0.9;
q=0.3;
p=1-q;

% initinalization

Price_of_Anarchy = zeros(1,length(Mu));
Price_of_S = zeros(1,length(Mu));
TR_1=zeros(1,length(Mu));
TR_2=zeros(1,length(Mu));
 Number_of_NE = zeros(1,length(Mu));
for kk = 1:1:length(Mu)
     
     mu=100+100.*(kk-1);
     
     TR = zeros(n_1,n_2); 
  TR_NE = zeros(n_1,n_2);
     
%% optimal solution    

     for ii= 1:1:(n_1+1)
       for jj=1:1:(n_2+1)
           
            m_1=ii-1;
            m_2=jj-1;
            
            T_1 = m_1*phi+(n_2-m_2)*p*phi;
            T_2 = m_2*phi+(n_1-m_1)*p*phi;
            TR(ii,jj) = (T_1*mu)/(T_1+mu)+(T_2*mu)/(T_2+mu);
        end
    end
    
TR_opt = max((max(TR)));
[opt_m1,opt_m2] = find(TR_opt==TR);
opt_m1=max(opt_m1);opt_m2=max(opt_m2);
% val with different q
val =  round((n_2 + n_1*p + (1-sqrt(p))*mu/phi)/(p + sqrt(p)));
val(val>n_1) = n_1;


%% Nash equilibrium 
  for ii= 1:1:(n_1+1)
     for jj=1:1:(n_2+1)
           
            m_1=ii-1;
            m_2=jj-1;

            t1= (q*mu/phi+m_2*(1+p^2)+(n_1+1)*p-n_2.*p^2)./(2*p);
            t2=(q*mu/phi+m_1*(1+p^2)+(n_2+1)*p-n_1.*p^2)./(2*p);
            
           % t1(t1==inf)=n1;
            %t2(t2==inf)=n2;
     % case 1
        if (m_1==n_1) && (m_2==n_2) && (m_1<=t1) && (m_2 <= t2)
            
            TR_NE(ii,jj) = mu.*(2 - mu/phi/(n_1+mu/phi) - mu/phi/(n_2+mu/phi));
            
        end   
        
     % case 2
        if m_2==n_2&& 0<m_1&& m_1<n_1 && (t1-1)<=m_1&& m_1<=t1
            
           TR_NE(ii,jj) = mu.*(2 - mu/phi/(m_1+mu/phi) - mu/phi/(n_2+(n_1-m_1)*p+mu/phi)); 

        end
        
     % case 5
        if m_1==0 && m_2==0 && m_1>=(t1-1)&& m_2>=(t2-1)
            
            TR_NE(ii,jj) = mu.*(2 - mu/phi/(n_2*p+mu/phi) - mu/phi/(n_1*p+(n_1-m_1)*p+mu/phi)); 
        end
       
     % case 6
        if 0<m_1&& m_1<n_1 && 0<m_2&& m_2<n_2 && m_1>=t1-1 && m_2>= (t2-1)
            
           TR_NE(ii,jj) = mu.*(2 - mu/phi/(m_2+(n_1-n_2)/2+(n_2-m_2)*p+mu/phi) - mu/phi/(m_2+((n_1+n_2)/2-m_2)*p + mu/phi)); 
        end 
        
     % case 8
        if 0<m_1&& m_1<n_1 && m_2==0 && t1>=m_1&& m_1>=(t1-1) && m_2>=(t2-1)
            
           TR_NE(ii,jj) = mu.*(2 - mu/phi/(m_2+(n_1-n_2)/2+(n_2-m_2)*p+mu/phi) - mu/phi/(m_2+((n_1+n_2)/2-m_2)*p + mu/phi)); 
        end 
        
     end
  end
  Number_of_NE(kk) =sum(sum(TR_NE~=0));
  %% Price_of_Anarchy
  %Price_of_S(kk) = TR_opt./ max(max(TR_NE));
  TR_NE(TR_NE==0)=inf;
  Price_of_Anarchy(kk) = TR_opt./ min(min(TR_NE));
  %% restore opt and ne
  TR_1(kk) = TR_opt;
  TR_2(kk)= min(min(TR_NE));
end

%% figure
figure(20)
hold on
plot(Mu,TR_1./TR_2,'*b-');
hold on
%plot(Mu,Price_of_S,'*k-');
xlabel('\mu');
ylabel('Price');
yyaxis right
plot(Mu,TR_1,'*k-');
hold on
plot(Mu,TR_2,'*k-');
legend('PoA','TR(NE)','TR(opt)');
grid on



