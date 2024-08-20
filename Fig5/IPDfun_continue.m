function dydt=IPDfun_continue(t,y,vaccine_varient,vaccine_time,vaccine_type,vaccine_amount,X,self_neu,parameter_set,pre_varient)

% ode of the model
%------parameters setting-------
% direct estimation

gam_Ag=parameter_set(4); %antigen degration rate
gam_Ab=parameter_set(5); %antibody degration
d_B=parameter_set(6); %naive B cell decay rate
d_M=parameter_set(7); %memory B cell decay rate

% primary immune parameters
s_pre=parameter_set(11);
s_F=parameter_set(12); %max GC maturation rate 
K=parameter_set(13); %Amount of antigen required to reach half the maximum maturation rate
p_B=parameter_set(14); %max naive antibody production rate
gam_ne=parameter_set(15); %antigen-antibody neutralization rate （baseline）

% memory immune parameters
p_M=parameter_set(16); %memory antibody production rate
k_B2M=parameter_set(17); % max GC B cell to memory B cell rate (germinal center)
max_M=parameter_set(21);

uni_varient=unique([pre_varient,vaccine_varient]);
uni_varient=uni_varient(uni_varient~=0);
n_varient=length(uni_varient); %number of unique vaccine type

%----load functions-----
Funs=EvaluationFunc(vaccine_varient,vaccine_type,vaccine_time,vaccine_amount,X,self_neu,parameter_set,n_varient);

gam_R=0;
k=0;

%R:y(1:n_type) Ag:y(n_type+1:2*n_type) Ab:y(2*n_type+1:3*n_type) F:y(3*n_type+1:4*n_type) M_off:y(4*n_type+1:5*n_type) M_on:y(5*n_type+1:6*n_type)
flux_dR=@(R) gam_R*R; %mRNA degration
flux_trans=@(R) k*R; %antigen translation
flux_neutral_Ag=@(Ag,Ab) gam_ne*(X*(Ab.*self_neu)).*Ag; %antibody antigen neutralizing
flux_neutral_Ab=@(Ag,Ab) 0; %抗原不导致抗体损失
flux_dAg=@(Ag) gam_Ag*Ag; %anigen degration
flux_B_Ab=@(F) p_B*F; %naive antibody production
flux_M_Ab=@(M) p_M*M; %memory antibody production
flux_dAb=@(Ab) gam_Ab*Ab; %antibody degration
flux_mature=@(Ag,F) s_F*Funs.saturation(Ag,K).*F.*(1-sum(F)); %GC B cell maturation
flux_F_decay=@(F) d_B*F; %GC B cell decay
flux_B2M=@(M_off,M_on,F) k_B2M*F*(1-sum(M_off+M_on)/max_M); %GC B cell differentiate into memory cells
% flux_M_prolif=@(Ag,M) s_M*((M.*X)./((M'*K)))*f(Ag,K);  %function
flux_M_decay=@(M) d_M*M; %memory B cell decay


flux_pre_mature=s_pre*Funs.a_pre(t,vaccine_time); %GC B cell pre mature

dR= -flux_dR(y(1:n_varient));

dAg= flux_trans(y(1:n_varient))-flux_neutral_Ag(y(n_varient+1:2*n_varient),y(2*n_varient+1:3*n_varient))-flux_dAg(y(n_varient+1:2*n_varient));

dAb= Funs.a(t,vaccine_time,y(n_varient+1:2*n_varient))*flux_B_Ab(y(3*n_varient+1:4*n_varient))+Funs.a_M(t,vaccine_time,y(n_varient+1:2*n_varient))*Funs.p_M_attune(t,vaccine_time).*flux_M_Ab(y(5*n_varient+1:6*n_varient))-flux_dAb(y(2*n_varient+1:3*n_varient))-flux_neutral_Ab(y(n_varient+1:2*n_varient),y(2*n_varient+1:3*n_varient));


% aaa=flux_pre_mature
% bbb=Funs.a(t,vaccine_time,y(n_varient+1:2*n_varient))*flux_mature(y(n_varient+1:2*n_varient),y(3*n_varient+1:4*n_varient))
% ccc=flux_F_decay(y(3*n_varient+1:4*n_varient))

dF=flux_pre_mature+Funs.a(t,vaccine_time,y(n_varient+1:2*n_varient))*flux_mature(y(n_varient+1:2*n_varient),y(3*n_varient+1:4*n_varient))-flux_F_decay(y(3*n_varient+1:4*n_varient));


dM_off=flux_B2M(y(4*n_varient+1:5*n_varient),y(5*n_varient+1:6*n_varient),y(3*n_varient+1:4*n_varient))-flux_M_decay(y(4*n_varient+1:5*n_varient));

dM_on=Funs.a_M(t,vaccine_time,y(n_varient+1:2*n_varient))*Funs.flux_M_prolif(y(n_varient+1:2*n_varient),y(4*n_varient+1:5*n_varient),y(5*n_varient+1:6*n_varient))-flux_M_decay(y(5*n_varient+1:6*n_varient));


dydt=[dR;dAg;dAb;dF;dM_off;dM_on];


end

