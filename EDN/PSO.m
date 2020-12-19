clc
clear
close all
format short g


%%%%%%%%%%% Uncompansated %%%%%%%%%%%%%%%%
start=0;
x=0;
[LSI1,LSI2,Vm,PTloss,QTloss,power_f_active]=pf_EDN(start,x);
disp(' ')
disp('======================================================')
disp('Results of EDN system without compansation')
disp(' ')
disp(['Total active loss is: ' num2str(PTloss) ' kW'])
Kp=168; % $/Kw.
T_cost=Kp*PTloss;
disp(['Total annual cost is: ' num2str(T_cost) ' $'])
[value_v,index_v]=sort(abs(Vm));
disp(['Minimum voltage is: ' num2str(value_v(1)) ', at bus ' num2str(index_v(1))])
disp(['Maximum voltage is: ' num2str(value_v(end-1)) ', at bus ' num2str(index_v(end-1))])
disp(' ')
pause(0.5)

figure(2);
plot(power_f_active)
%%%%%%%%%%% Compansation %%%%%%%%%%%%%%%%%%
%% %%%%%%%%% Case 1 %%%%%%%%%%%%%%%%%%%%%%%%
start=1;
NC_max=4;
%% parameters setting
a=numel(LSI1(:,1))/2;
a=round(a);
lb=[LSI1(1:a,2)' 0*ones(1,a)]; % lower bound
ub=[LSI1(1:a,2)' 1200*ones(1,a)];  % upper bound
nvar=2*a; % number of variable


NP=500;              % number particle
T=100;               % max of iteration

W=1;
C1=2;
C2=2;

alpha=0.05;

%% initialization
tic
empty.pos=[];
empty.cost=[];
empty.velocity=[];

particle=repmat(empty,NP,1);

for i=1:NP
particle(i).pos=lb+rand(1,nvar).*(ub-lb);
[particle(i).cost]=fitness(NC_max,particle(i).pos,start);
particle(i).velocity=0;
end

bparticle=particle;

[value,index]=min([particle.cost]);

gparticle=particle(index);

%% main loop

best=zeros(T,1);
AVR=zeros(T,1);

for t=1:T

     for i=1:NP
         
          particle(i).velocity=W*particle(i).velocity...
                              +C1*rand(1,nvar).*(bparticle(i).pos-particle(i).pos)...
                              +C2*rand(1,nvar).*(gparticle.pos-particle(i).pos);
          
         particle(i).pos=particle(i).pos+particle(i).velocity;
         
         particle(i).pos=min(particle(i).pos,ub);
         particle(i).pos=max(particle(i).pos,lb);
          
         
         [particle(i).cost]=fitness(NC_max,particle(i).pos,start);
         
         
         if particle(i).cost<bparticle(i).cost
             bparticle(i)=particle(i);
             
             if bparticle(i).cost<gparticle.cost
                 gparticle=bparticle(i);
             end
         end
         
         

     end
     
     
     
 W=W*(1-alpha);
 
 best(t)=gparticle.cost;
 AVR(t)=mean([particle.cost]);
 
 disp([ ' t = ' num2str(t)   ' BEST = '  num2str(best(t))]);
 

 
end

%% results
disp('====================================================')
disp([' Buses are   =  '  num2str(gparticle.pos(1:a))])
disp([' sizes are   =  '  num2str(gparticle.pos(a+1:end))])

j=0;
for i=1:a
    if gparticle.pos(a+i)<150
    gparticle.pos(a+i)=0;
    else
    j=j+1;
    Buses(j)=gparticle.pos(i);
    Sizes(j)=gparticle.pos(a+i);
    end
end
    disp(' ')
    disp('          Buses    Size(kVAr) ')
    disp([Buses' Sizes'])   
    disp(' ')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LSI1,LSI2,Vm,PTloss,QTloss,power_f_active]=pf_EDN(start,gparticle.pos);

figure(1)
hold on
plot(abs(Vm),'-r^')

figure(2);
hold on
plot(power_f_active,'-r^')

disp(' ')
disp('======================================================')
disp('Results of EDN system with compansation, Case 1')
disp(' ')
disp(['Total active loss is: ' num2str(PTloss) ' kW'])
Kp=168; % $/Kw
Kc=5; % $/Kvar
life_exp=10; %life expectancy
T_cost=Kp*PTloss;
disp(['Total annual cost is: ' num2str(T_cost) ' $'])
[value_v,index_v]=sort(abs(Vm));
disp(['Minimum voltage is: ' num2str(value_v(1)) ', at bus ' num2str(index_v(1))])
disp(['Maximum voltage is: ' num2str(value_v(end-1)) ', at bus ' num2str(index_v(end-1))])
Total_cap_cost=Kc*sum(gparticle.pos(a+1:end))/life_exp;
disp(['Total Capacitor cost is: ' num2str(Total_cap_cost) ' $'])
disp(' ')

%% %%%%%%%%%%%% Case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%

start=2;
NC_max=4;
%% parameters setting
a=numel(LSI1(:,1))/2;
a=round(a);
disp(' ')
disp('=======================================')
disp('           Case 2           ')
lb=[LSI1(1:a,2)' 0*ones(1,a)]; % lower bound
ub=[LSI1(1:a,2)' 1*ones(1,a)];  % upper bound

nvar=2*a; % number of variable
NP=500;              % number particle
T=100;               % max of iteration
W=1;
C1=2;
C2=2;
alpha=0.05;

%% initialization
tic
empty.pos=[];
empty.cost=[];
empty.velocity=[];

particle=repmat(empty,NP,1);

for i=1:NP
particle(i).pos=lb+rand(1,nvar).*(ub-lb);
[particle(i).cost]=fitness(NC_max,particle(i).pos,start);
particle(i).velocity=0;
end

bparticle=particle;

[value,index]=min([particle.cost]);

gparticle=particle(index);

%% main loop

best=zeros(T,1);
AVR=zeros(T,1);

for t=1:T

     for i=1:NP
         
          particle(i).velocity=W*particle(i).velocity...
                              +C1*rand(1,nvar).*(bparticle(i).pos-particle(i).pos)...
                              +C2*rand(1,nvar).*(gparticle.pos-particle(i).pos);
          
         particle(i).pos=particle(i).pos+particle(i).velocity;
         
         particle(i).pos=min(particle(i).pos,ub);
         particle(i).pos=max(particle(i).pos,lb);
          
         
         [particle(i).cost]=fitness(NC_max,particle(i).pos,start);
         
         
         if particle(i).cost<bparticle(i).cost
             bparticle(i)=particle(i);
             
             if bparticle(i).cost<gparticle.cost
                 gparticle=bparticle(i);
             end
         end
         
         

     end
     
     
     
 W=W*(1-alpha);
 
 best(t)=gparticle.cost;
 AVR(t)=mean([particle.cost]);
 
 disp([ ' t = ' num2str(t)   ' BEST = '  num2str(best(t))]);
 

 
end

%% results
disp('====================================================')
disp([' Buses are   =  '  num2str(gparticle.pos(1:a))])
disp([' sizes are   =  '  num2str(gparticle.pos(a+1:end))])

x=gparticle.pos;
for i=1:numel(x)/2
if x(i+(numel(x)/2))>0 && x(i+(numel(x)/2))<=0.125
            x(i+(numel(x)/2))=150;
        elseif x(i+(numel(x)/2))>0.125 && x(i+(numel(x)/2))<=0.25
            x(i+(numel(x)/2))=300;
        elseif x(i+(numel(x)/2))>0.25 && x(i+(numel(x)/2))<=0.375
            x(i+(numel(x)/2))=450;
        elseif x(i+(numel(x)/2))>0.375 && x(i+(numel(x)/2))<=0.5
            x(i+(numel(x)/2))=600;
        elseif x(i+(numel(x)/2))>0.5 && x(i+(numel(x)/2))<=0.625
            x(i+(numel(x)/2))=750;
        elseif x(i+(numel(x)/2))>0.625 && x(i+(numel(x)/2))<=0.75
            x(i+(numel(x)/2))=900;
        elseif x(i+(numel(x)/2))>0.75 && x(i+(numel(x)/2))<=0.875
            x(i+(numel(x)/2))=1050;
        elseif x(i+(numel(x)/2))>0.875
            x(i+(numel(x)/2))=1200;
end
end

j=0;
for i=1:a
    if x(a+i)~=0
    j=j+1;
    Buses(j)=x(i);
    Sizes(j)=x(a+i);
    end
end
    disp(' ')
    disp('          Buses    Size(kVAr) ')
    disp([Buses' Sizes'])   
    disp(' ')
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[LSI1,LSI2,Vm,PTloss,QTloss,power_f_active]=pf_EDN(start,gparticle.pos);

figure(1)
hold on
plot(abs(Vm),'-k*')
legend('Without Compansation','case1','case2')
xlim([1 30])

figure(2);
hold on
plot(power_f_active,'-k*')
legend('Without Compansation','case1','case2')
xlim([1 30])

disp(' ')
disp('======================================================')
disp('Results of EDN system with compansation, Case 2')
disp(' ')
disp(['Total active loss is: ' num2str(PTloss) ' kW'])
Kp=168; % $/Kw
Kc=5; % $/Kvar
life_exp=10; %life expectancy
T_cost=Kp*PTloss;
disp(['Total annual cost is: ' num2str(T_cost) ' $'])
[value_v,index_v]=sort(abs(Vm));
disp(['Minimum voltage is: ' num2str(value_v(1)) ', at bus ' num2str(index_v(1))])
disp(['Maximum voltage is: ' num2str(value_v(end-1)) ', at bus ' num2str(index_v(end-1))])
Total_cap_cost=Kc*sum(x(a+1:end))/life_exp;
disp(['Total Capacitor cost is: ' num2str(Total_cap_cost) ' $'])
disp(' ')


save PSO
