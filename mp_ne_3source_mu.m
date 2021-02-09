%%  multi-cource,pure strategy, general case


%clc;
%clear;
%close all;

%% parameter seeting

%traffic generation rate
phi = 1;

%traffic average service time
Mu = 0:1:10;

% number of sources
m = 3;

% number of players at each node
%N_m=sort(randperm(3)*4,'descend');
N_m=[12,8,4];
n = sum(N_m);
n_1 = N_m(1);
n_2 = N_m(2);
n_3 = N_m(3);

% IP loss probability
q=0.5;
p = 1-q;
TR_NE1=zeros(1,length(Mu));
for kk=1:length(Mu)
mu=kk-1;
% initinalization
number_of_NE=0;
A=zeros(6,200);

for u1 = 0:n_1  
    for u2 = 0:n_2
        for u3 = 0:n_3
            for v1 = 0:n-u1-u2-u3
                for v2 = 0:n-u1-u2-u3-v1
                    for  v3 = n-u1-u2-u3-v1-v2
                        for x_12 = 0:n_1-u1      
                            for x_13 = 0:n_1-u1-x_12
                                for x_21 = 0:n_2-u2
                                    for x_23 = 0:n_2-u2-x_21
                                        for x_31 = 0:n_3-u3
                                            for x_32 = 0:n_3-u3-x_31
                                                
%% if truth, continue                                                
if x_12+x_32==v2 && x_21+x_31==v1 && x_13+x_23==v3
                                                  
                                                    
                                                    %% main part
                                                    % u_star,v_star
                                                    au1 = u1+p*v1; au2 = u2+p*v2; au3 = u3+p*v3;
                                                    if au1==min(au1,min(au2,au3))
                                                        u_star = u1; v_star=v1;
                                                    elseif au2 == min(au1,min(au2,au3))
                                                        u_star = u2; v_star=v2;
                                                    elseif au3 == min(au1,min(au2,au3))
                                                        u_star = u3; v_star=v3;
                                                    end
                                                    % initialization
                                                    c_19_1=[];c_19_2=[];c_19_3=[];
                                                    c_20_1=[];c_20_2=[];c_20_3=[];c_20_4=[];c_20_5=[];c_20_6=[];
                                                    c_21_1=[];c_21_2=[];c_21_3=[];c_21_4=[];c_21_5=[];c_21_6=[];
                                                    
                                                    % check condition (19)
                                                    B_19= u_star+p*(v_star+1);
                                                    if u1>0
                                                        c_19_1 = p*(u1+p.*v1)-q*mu/phi-B_19;
                                                    end
                                                    
                                                    if u2>0
                                                        c_19_2 = p*(u2+p.*v2)-q*mu/phi-B_19;
                                                    end
                                                    
                                                    if u3>0
                                                        c_19_3 = p*(u3+p.*v3)-q*mu/phi-B_19;
                                                    end
                                                    
                                                    % check condition(20)
                                                    if x_12>0          
                                                        c_20_1 = u2+v2*p-p*(u1+1+p*v1)+q*mu/phi;
                                                        c_21_1 = u2+v2*p-u_star-p*(v_star+1);
                                                    end
                                                    if x_13>0
                                                        c_20_2= u3+v3*p-p*(u1+1+p*v1)+q*mu/phi;
                                                        c_21_2 = u3+v3*p-u_star-p*(v_star+1);
                                                    end
                                                    if x_21>0
                                                        c_20_3= u1+v1*p-p*(u2+1+p*v2)+q*mu/phi;
                                                        c_21_3 = u1+v1*p-u_star-p*(v_star+1);
                                                    end
                                                    if x_23>0
                                                        c_20_4= u3+v3*p-p*(u2+1+p*v2)+q*mu/phi;
                                                        c_21_4 = u3+v3*p-u_star-p*(v_star+1);
                                                    end
                                                    if x_31>0
                                                        c_20_5= u1+v1*p-p*(u3+1+p*v3)+q*mu/phi;
                                                        c_21_5 = u1+v1*p-u_star-p*(v_star+1);
                                                    end
                                                    if x_32>0
                                                        c_20_6= u2+v2*p-p*(u3+1+p*v3)+q*mu/phi;
                                                        c_21_6 = u2+v2*p-u_star-p*(v_star+1);
                                                    end
                                                    
                                                    %check logic
                                                    
                                                    condition = [c_19_1 c_19_2 c_19_3 c_20_1 c_20_2 c_20_3...
                                                                c_20_4 c_20_5 c_20_6 c_21_1 c_21_2 c_21_3 c_21_4 c_21_5 c_21_6];
                                                    logics = sum(find(condition>0));
                                                    
                                                    output=[];
                                                    if logics==0
                                                        number_of_NE=number_of_NE+1;
                                                        output = [u1 u2 u3 v1 v2 v3];
                                                        A(:,number_of_NE)=output;
                                                       
                                                     end
                                                    
                                                    %% main part end
end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
 A=A(:,1:number_of_NE);
 TR=zeros(1,number_of_NE);
for ii=1:number_of_NE
TR(ii)= mu*(m-(mu./(A(1,ii)*phi+A(4,ii).*p*phi+mu)+mu./(A(2,ii)*phi+A(5,ii).*p*phi+mu)+mu./(A(3,ii)*phi+A(6,ii).*p*phi+mu)));
end
TR_NE1(kk)= min(TR);
end
%-----------------------------------------------
%TR(ii)= mu*(m-(mu./(A(1)*phi+A(4).*p*phi+mu)+mu./(A(2)*phi+A(5).*p*phi+mu)+mu./(A(3)*phi+A(6).*p*phi+mu)));
%for ii=1:27
%TR(ii)= mu*(m-(mu./(A(1,ii)*phi+A(4,ii).*p*phi+mu)+mu./(A(2,ii)*phi+A(5,ii).*p*phi+mu)+mu./(A(3,ii)*phi+A(6,ii).*p*phi+mu)));
%end