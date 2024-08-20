function [t,Rt,Agt,Abt,Ft,M_offt,M_ont]=simIPD_continue(vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,X,self_neu,parameter_set,gap_sample,pre_t,pre_Rt,pre_Agt,pre_Abt,pre_Ft,pre_M_offt,pre_M_ont,pre_varient)
%------parameters setting-------
% direct estimation

Pug=parameter_set(1); % antigen amount per ug
PRug=parameter_set(2);% antigen amount per ug(mRNA vaccine)
Vug=parameter_set(3); % antigen amount per ug(灭活病毒疫苗)

n_vaccine=length(vaccine_varient); %number of vaccine dose

uni_varient=unique([pre_varient,vaccine_varient]);
uni_varient=uni_varient(uni_varient~=0);
n_varient=length(uni_varient); %number of unique vaccine type

% if sum(vaccine_type==1)
%     self_neu(1)=parameter_set(20); %WT antibody self neutrlazing relative to BA.5
% end

for vi=(1:n_varient)
    v=uni_varient(vi);
    vaccine_varient(vaccine_varient==v)=vi;
end

p_Rt=zeros(n_varient,length(pre_t));
p_Rt(pre_varient,:)=pre_Rt;

p_Agt=zeros(n_varient,length(pre_t));
p_Agt(pre_varient,:)=pre_Agt;

p_Abt=zeros(n_varient,length(pre_t));
p_Abt(pre_varient,:)=pre_Abt;

p_Ft=zeros(n_varient,length(pre_t));
p_Ft(pre_varient,:)=pre_Ft;

p_M_offt=zeros(n_varient,length(pre_t));
p_M_offt(pre_varient,:)=pre_M_offt;

p_M_ont=zeros(n_varient,length(pre_t));
p_M_ont(pre_varient,:)=pre_M_ont;



options=odeset('RelTol',0.000001); %dde numerical solution option

Rt=[];
Agt=[];
Abt=[];
Ft=[];
M_offt=[];
M_ont=[];
t=[];
%------run simulation-------
for dose=(1:n_vaccine)
    if dose==n_vaccine
        tspan=[vaccine_time(dose),vaccine_time(dose)+gap_sample]; %the last dose
    else
        tspan=[vaccine_time(dose),vaccine_time(dose+1)]; %simulation betwwen 2 dose
    end
    
    if dose==1
        y0=zeros(6*n_varient,1);
        y0(pre_varient)=pre_Rt(:,end);
        y0(pre_varient+n_varient)=pre_Agt(:,end);
        y0(pre_varient+2*n_varient)=pre_Abt(:,end);
        y0(pre_varient+3*n_varient)=pre_Ft(:,end);
        y0(pre_varient+4*n_varient)=pre_M_offt(:,end);
        y0(pre_varient+5*n_varient)=pre_M_ont(:,end);
    else
        y0=sol.y(:,end); %last solution as the y0

    end
    injection=vaccine_varient(:,dose);
    injection=injection(injection~=0);
    valent=length(injection); %the value of vaccine
    for injection_i=injection
        if vaccine_type(dose)==1 % the vaccine type is antigen 
            if injection_i==1 %the vaccine is CoronaVac
                P0=Vug*vaccine_amount(dose)/valent; %无论单价，多价，每一针的总剂量一致。
            else
                P0=Pug*vaccine_amount(dose)/valent;
            end
            y0(injection_i+n_varient)=y0(injection_i+n_varient)+P0; %antiegn injection
        elseif vaccine_type(dose)==2 % the vaccine type is mRNA
            P0=PRug*vaccine_amount(dose)/valent;
            y0(injection_i+n_varient)=y0(injection_i+n_varient)+P0; %mRNA injection
        end
    end
    
    y0(5*n_varient+1:6*n_varient)=y0(5*n_varient+1:6*n_varient)+y0(4*n_varient+1:5*n_varient);
    y0(4*n_varient+1:5*n_varient)=0;
    sol = ode23(@(t,y) IPDfun_continue(t,y,vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,X,self_neu,parameter_set,pre_varient),tspan,y0,options);
    
    %-------get output-------
    Rt=[Rt,sol.y(1:n_varient,1:end-1)];
    Agt=[Agt,sol.y(n_varient+1:2*n_varient,1:end-1)];
    Abt=[Abt,sol.y(2*n_varient+1:3*n_varient,1:end-1)];
    Ft=[Ft,sol.y(3*n_varient+1:4*n_varient,1:end-1)];
    M_offt=[M_offt,sol.y(4*n_varient+1:5*n_varient,1:end-1)];
    M_ont=[M_ont,sol.y(5*n_varient+1:6*n_varient,1:end-1)];
    t=[t,sol.x(1:end-1)];
end
    Rt=[p_Rt(:,1:end-1),Rt];
    Agt=[p_Agt(:,1:end-1),Agt];
    Abt=[p_Abt(:,1:end-1),Abt];
    Ft=[p_Ft(:,1:end-1),Ft];
    M_offt=[p_M_offt(:,1:end-1),M_offt];
    M_ont=[p_M_ont(:,1:end-1), M_ont];
    t=[pre_t(:,1:end-1),t];
end