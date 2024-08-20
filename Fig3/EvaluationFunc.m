function Funcollect=EvaluationFunc(vaccine_varient,vaccine_type,vaccine_time,vaccine_amount,X,self_neu,parameter_set)
Funcollect.flux_M_prolif=@(Ag,M_off,M_on) flux_M_prolif(Ag,M_off,M_on,X,self_neu,parameter_set);
Funcollect.saturation=@saturation;
Funcollect.a=@(t,vaccine_time,Ag) a(t,vaccine_time,Ag,parameter_set);
Funcollect.a_M=@(t,vaccine_time,Ag) a_M(t,vaccine_time,Ag,parameter_set);
Funcollect.a_pre=@(t,vaccine_time) a_pre(t,vaccine_varient,vaccine_time,parameter_set);
Funcollect.ddehist=@(t) ddehist(t,vaccine_varient,vaccine_type,vaccine_time,vaccine_amount,parameter_set);
Funcollect.p_M_attune=@(t,vaccine_time) p_M_attune(t,vaccine_varient,vaccine_time,X,parameter_set);

end
% 自定义函数集合
function flux=flux_M_prolif(Ag,M_off,M_on,X,self_neu,parameter_set) %affinity dependent memory B cell proliferation
s_M=parameter_set(18);
K=parameter_set(13);
max_M=parameter_set(21); %max Memory B cell level
M_n=parameter_set(19);

if sum(M_on)
    flux=s_M*((((M_on.*self_neu).*X')./(((M_on.*self_neu)'*X'+M_n)))*saturation(Ag,K)).*M_on*(1-sum(M_on+M_off)/max_M);
else
    flux=0;
end
end

function DC=saturation(Ag,K)
DC=Ag./(Ag+K);
end

function active=a(t,vaccine_time,Ag,parameter_set) 
t_start=parameter_set(10);
t_max=parameter_set(9);
t_min=parameter_set(8);

judge1=(vaccine_time+t_start<t) & (t<=vaccine_time+t_min);%immune response time span

judge2=(vaccine_time+t_min<t) & (t<=vaccine_time+t_max);

active=(sum(judge1)>0)||((sum(judge2)>0)&&(sum(Ag)>1));

end

function active=a_M(t,vaccine_time,Ag,parameter_set) 
t_max=parameter_set(9);
t_min=parameter_set(8);

judge1=(vaccine_time<t) & (t<=vaccine_time+t_min);%immune response time span

judge2=(vaccine_time+t_min<t) & (t<=vaccine_time+t_max);

active=(sum(judge1)>0)||((sum(judge2)>0)&&(sum(Ag)>1));

end

function active=a_pre(t,vaccine_varient,vaccine_time,parameter_set) 
t_start=parameter_set(10);
uni_varient=unique(vaccine_varient);
uni_varient=uni_varient(uni_varient~=0);
n_varient=length(uni_varient); %number of unique vaccine type
active=zeros(n_varient,1);
judge=(vaccine_time<t) & (t<vaccine_time+t_start);%naive B cell pre maturation time span
if sum(judge)>0
    ind=judge==1; %find which dose
    type=vaccine_varient(:,ind); %get the vaccine type
    type=type(type~=0);
    active(type)=1; % active the corresponding type GC B cell pre maturation
end
    
end


function attune=p_M_attune(t,vaccine_varient,vaccine_time,X,parameter_set)
Ab_n=parameter_set(22);
judge=(t>=vaccine_time)&(t<[vaccine_time(2:end),10000]);%生成当前在第几次注射后
injection=vaccine_varient(:,judge);
injection=injection(injection~=0);
cross_neu=X(injection,:);

cross_neu=max(cross_neu,[],1)';%记忆B细胞与当前抗体最大的交叉免疫
%cross_neu=max(cross_neu)';
attune=cross_neu./(cross_neu+Ab_n);

end





function s = ddehist(t,vaccine_varient,vaccine_type,vaccine_time,vaccine_amount,parameter_set)
% history function for dde simulation

Pug=parameter_set(1);
PRug=parameter_set(2);
Vug=parameter_set(3); % antigen amount per ug
n_varient=length(unique(vaccine_varient));
s = zeros(5*n_varient,1);
if t==vaccine_time(1) 
    injection=vaccine_varient(1);
    if vaccine_type(1)==1 % the vaccine type is antigen 
        if vaccine_varient(1)==1 %the vaccine is CoronaVac
            P0=Vug*vaccine_amount(1);
        else
            P0=Pug*vaccine_amount(1);
        end
        s(injection+n_varient)=P0; %antiegn injection
    elseif vaccine_type(1)==2 % the vaccine type is mRNA
        P0=PRug*vaccine_amount(1);
        s(injection+n_varient)=P0; %mRNA injection
    end    
end
end