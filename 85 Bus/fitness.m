function Z=fitness(NC_max,x,start)

[LSI1,LSI2,V_sys,PTloss,QTloss,power_f_active]=pf_85(start,x);

Kp=168; % $/Kw.
Kc=5; % $/Kvar
life_exp=10; %life expectancy
Z=Kp*PTloss+Kc*sum(x(1+(numel(x)/2):end));

%%%%% voltage constraint
if max(abs(V_sys))>1.1 || min(abs(V_sys))<0.9
Z=sum(abs(V_sys))*Z;
end

%%%%% 
a=x(round(numel(x)/2)+1:end)>0;
if sum(a)>NC_max 
    Z=sum(a)*Z;
end

end
