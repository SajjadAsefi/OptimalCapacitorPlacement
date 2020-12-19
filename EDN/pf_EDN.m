function [LSI1,LSI2,Vm,PTloss,QTloss,power_f_active]=pf_EDN(start,x)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Network_data_EDN;

S_base=27221;   %(kVA)
V_base=11; %(kV)
Z_base=1000*(V_base^2)/(S_base);
I_base=S_base/(sqrt(3)*V_base);
if start==1
for i=1:numel(x)/2
    bus_data(x(i),3)=bus_data(x(i),3)-x(i+(numel(x)/2));   
end
elseif start==2
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
for i=1:numel(x)/2
 bus_data(x(i),3)=bus_data(x(i),3)-x(i+(numel(x)/2)); 
end

end
demanded_P=bus_data(:,2)/S_base;
demanded_Q=bus_data(:,3)/S_base;
R=line_data(:,3)/Z_base;
X=line_data(:,4)/Z_base;
nbus=length(bus_data);
nline=length(line_data);
%% calculating configuration matrix of network

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% uncomment in case you want to view graphically the system
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from=line_data(:,1);
% to=line_data(:,2);

% uncomment in case you want to draw the system figure
% a=from;
% b=to;
% w=line_data(:,3);
% u=max(max(a),max(b));
% DG = sparse(a,b,w,u,u); 
% pathMAT=graphallshortestpaths(DG);
% for i=1:nbus
%     qq(i)={num2str(i)};
% end
% if start == 0
%   view(biograph(DG,qq,'ShowWeights','off'))
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BIBC=zeros(nbus-1,nline);
for i=1:nline
    BIBC(:,line_data(i,2))=BIBC(:,line_data(i,1));
    BIBC(line_data(i,5),line_data(i,2))=1;
end
BIBC(:,1)=[]; %BIBC Matrix

%% Initilize voltage 
Vm=ones(1,nbus);   %voltage magnitude(pu)

%% start iteration
delta=1;eps=0.00001;iter=0;
MAXiter=1000;
while delta > eps
    iter=iter+1;
    if iter>MAXiter
        break;
    end
    
    for k=1:nbus
        Ibus(k,1)=(demanded_P(k)-sqrt(-1).*demanded_Q(k))./(conj(Vm(k)));
    end
    Inode=BIBC*Ibus(2:end); % It gives out branch current

    V(1)=1;
    for k=1:length(line_data)
       V(line_data(k,2))=Vm(line_data(k,1))-(R(k,1)+sqrt(-1)*X(k,1))*(Inode(line_data(k,5)));
    end
    delta=max(abs(V-Vm));
    Vm=V;
end
%% results 
if start==0
    figure(1)
    plot(abs(Vm),'--*b')
    ylabel('Voltage Magnitude (p.u.)')
    xlabel('Bus Number')
    title('Vlotage profile')
end
%% loss calculation
    IL=Inode;
for k=1:length(line_data)
%     IL(k)=(V(line_data(k,1))-V(line_data(k,2)))/(R(k)+sqrt(-1).*X(k));
    S_line(k)=V(line_data(k,1))*conj(IL(k));
    PLoss(k)=R(k)*(abs(IL(k))^2)*S_base;
    QLoss(k)=X(k)*(abs(IL(k))^2)*S_base;
end

power_f_active=real(S_line)*S_base;
power_f_reactive=imag(S_line)*S_base;

PTloss=sum(PLoss); % total active or reactive loss
QTloss=sum(QLoss);% total reactive loss

%%%%%%%%%%%%%%% loss sensitivity index %%%%%%%%%%%%
LSI1=-2*((demanded_P(2:end)+PLoss'/S_base).^2+(demanded_Q(2:end)+QLoss'/S_base).^2).*(R(k,1).*Z_base)./(abs(Vm(1,2:end)).^3)';
% LSI1=-2*(demanded_P(2:end).^2+demanded_Q(2:end).^2).*(R(k,1).*Z_base)./(abs(Vm(1,2:end)).^3)';
[Value1,Index1]=sort(LSI1);
Index1=Index1+1;
LSI1=[Value1,Index1];
    
LSI2=2*((demanded_Q(2:end)+QLoss'/S_base).^2).*(R(k,1).*Z_base)./(abs(Vm(1,2:end)).^2)';
[Value2,Index2]=sort(LSI2,'descend');
Index2=Index2+1;
LSI2=[Value2,Index2];