clc,clear,close all

%-----input batch vaccination-------
load('batch_information_human.mat')
batch_sample(:,end)=[];
batch_sample(2:3,end)=0;

load('X_yunda.mat')
load('self_neu_yunda.mat');
self_neu=self_neu';
load('vaccine_correspond_yunda.mat')
load('X_pse_yunda.mat')
load('titers_human_all.mat')
titers{1}(:,end)=[];
titers{2}(:,end)=[];
titers{3}(:,end)=[];

titers{1}(end,:)=[];
titers{2}(end,:)=[];
titers{3}(end,:)=[];

%----titer_correction--------
load('mix_parameters.mat')
X(5,1)=X(5,1)*sqrt(parameter_set(25));
X(6,1)=X(6,1)*parameter_set(28);
X(6,2)=X(6,2)*parameter_set(30);
X(6,3)=X(6,3)*parameter_set(31);
X(6,4)=X(6,4)*parameter_set(32);

X_pse(7:9,1)=X_pse(7:9,1).*parameter_set(25:27);
X_pse(11:12,1)=X_pse(11:12,1).*parameter_set(28:29);
X_pse(11,2:4)=X_pse(11,2:4).*parameter_set(30:32)';


pse={'WT','Alpha','Beta','Delta','BA.1','BA.2','BA.5','BF.7','BQ.1.1','XBB.1','XBB.1.5','JN.1'};
sub_pse=[1,2,3,7,8,9,11];

X_pse=X_pse(sub_pse,:);

%-----load parameters-------
load('human_parameter.mat')
parameter_set_new=parameter_set(1:22);
%parameter_set_new([2,20])=parameter_set(23:24);


% parameter_set_new(6)=parameter_set_new(6)*100;
%--------run simulation ----------

[batch_t,batch_Rt,batch_Agt,batch_Abt,batch_Ft,batch_Mofft,batch_Mont]=batch_simIPD(batch_vaccine_varient,batch_vaccine_time,batch_vaccine_type,batch_vaccine_amount,X,self_neu,parameter_set_new,vaccine_correspond);

[batch_titer_lgFC,batch_titer_est_lgFC,titer_group,titer_time,titer_sample]=batch_getFC_human(batch_vaccine_varient,batch_vaccine_time,batch_t,batch_Abt,batch_sample,titers,X_pse,self_neu,vaccine_correspond);
weight=[1,1,1];
%batch_titer_lgFC=batch_titer_lgFC(1:end-1,:);

[loss_weight,n_inverse]=batch_getloss_weight_human(batch_titer_lgFC,batch_titer_est_lgFC,weight);


%%
clc,close all

X = imread('fig4a.png');
w_X=size(X,2);
l_X=size(X,1);

n_pse=size(X_pse,1);
figure(1)
set(gcf,'Position',[50,50,800,950],'Color','w')
f_pos=get(gcf,'Position');
w_f=f_pos(3);
l_f=f_pos(4);

l_jj_Rsq=0.08;
l_b_up=0.015;
l_b=0.05;

w_b_r=0.03;
w_b=0.07;
w_b_top=0.05;


w_hum=1-2*w_b_top;
l_hum=w_hum*w_f*l_X/(l_f*w_X);
l_hum_Rsq=0.03;

l_Rsq=(1-l_b_up-l_b-l_jj_Rsq-l_hum-l_hum_Rsq)/2;

w_Rsq=l_Rsq*f_pos(4)/f_pos(3);

w_jj_leg_Rsq=0.1;

l_jj_pse=0.03;
l_comb=0.175;

l_pse=(1-l_comb-l_b_up-l_b-l_hum-l_hum_Rsq-(n_pse-1)*l_jj_pse-l_jj_pse)/n_pse;
w_leg=0.05;
w_jj_leg=0.02;
w_jj_leg_pse=0;

l_b_leg=0.06;
w_left=1-w_b_r-w_Rsq-w_jj_leg_Rsq-w_leg-w_jj_leg-w_b;
l_leg_pse=0.05;

%--------human_vaccination-------

ax_hum=axes;
set(ax_hum,'Position',[w_b_top,1-l_b_up-l_hum,w_hum,l_hum])


image(X(:,:,1:3))
axis off

%----------------%
n_titer=size(batch_titer_lgFC,2);

labels=pse(sub_pse);

order_index=(1:n_titer);
alpha=0.35;

c_exp=[255 165 0]/255;
% c_est=[107 62 151]/255;
c_est=[0 0 139]/255;

for i=(1:size(X_pse,1))
    ax_pse(i)=axes;
    set(ax_pse(i),'Position',[w_b,1-l_b_up-l_hum-l_hum_Rsq-l_pse-(i-1)*(l_pse+l_jj_pse),w_left,l_pse])

   
    hold on

    
    scatter((1:n_titer),batch_titer_lgFC(i,order_index),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_exp,'MarkerEdgeColor',c_exp,'LineWidth',1)
    scatter((1:n_titer),batch_titer_est_lgFC(i,order_index),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_est,'MarkerEdgeColor',c_est,'LineWidth',1)
    plot((1:n_titer),batch_titer_lgFC(i,order_index),'-','Color',c_exp)
    plot((1:n_titer),batch_titer_est_lgFC(i,order_index),'-','Color',c_est)
    
    
    xlabel(labels(i))

    
    min_FC=min([batch_titer_lgFC(i,:),batch_titer_est_lgFC(i,:)],[],2,'includenan');
    max_FC=max([batch_titer_lgFC(i,:),batch_titer_est_lgFC(i,:)],[],2,'includenan');

    min_FC=min(min_FC);
    max_FC=max(max_FC);
    
    min_FC=min(min_FC,-1.8);
    max_FC=max(max_FC,2.5);
    
    lgFC_width=[min_FC-1,max_FC+1];
    ylim(lgFC_width)
    set(gca,'YTick',(round(min_FC):2:round(max_FC)))
    
    
    ax_pse(i).FontSize=5;
    ax_pse(i).XLabel.FontSize=8;
    
    if i==4
        ylabel('log_2(titer FC)','FontSize',10)
    end
    xlim([0.5,n_titer+0.5])
    set(gca,'XTick',(1:n_titer))
    set(gca,'TickDir','out')
end

ax_comb=axes;
set(ax_comb,'Position',[w_b,l_b,w_left,l_comb])
hold on
varient_color=[127,127,127;247,187,127;16,131,88;194,127,159;85,160,251;198,55,53]/255;

num=0;
y_space=6;

for i=order_index
    num=num+1;
    group=titer_group(i);
    sample=titer_sample(i);
    time=titer_time(i);%采样时间
    
    vaccine_varient_i= batch_vaccine_varient(group,1:sample);
    vaccine_varinet_i=vaccine_correspond(vaccine_varient_i,:)';
    
    
    vaccine_time_i= batch_vaccine_time(group,1:sample);
    vaccine_type_i= batch_vaccine_type(group,1:sample);
    vaccine_amount_i= batch_vaccine_amount(group,1:sample);
    vaccine_interval=vaccine_time_i(2:end)-vaccine_time_i(1:end-1);
    
    interval_str=cell(1,sample);
    %size_i=scatter_sizes(1).*(vaccine_amount_i==1)+scatter_sizes(2).*(vaccine_amount_i==3)+scatter_sizes(3).*(vaccine_amount_i==10);
    
    %type_pr=(vaccine_type_i==1);
    %type_mrna=(vaccine_type_i==2);
    
    for j=(1:sample-1)
        interval_str{j}=num2str(vaccine_interval(j));
    end
    interval_str{sample}=num2str(time-vaccine_time_i(end));%最后一次注射与sample之间的间隔
    
    %XX=num*ones(1,sample);
    %YY=(110:-20:110-20*sample);
    
    for k=(1:sample)
        injection=vaccine_varinet_i(:,k);
        injection=injection(injection~=0);
        XX=num*ones(size(injection));
        valent=length(injection);
        y_center=130-20*k;
        YY=(y_center-(valent-1)*y_space/2:y_space:y_center+(valent-1)*y_space/2);
        CC=varient_color(injection,:);
        if injection==1
            scatter(XX,YY,30,CC,'filled')
        else
            scatter(XX,YY,80,CC,'s','filled')
        end
        if valent==2
            y_white=7;
            YY=y_center+[-y_white,y_white];
            scatter(XX,YY,80,'s','filled','MarkerFaceColor','w','MarkerEdgeColor','none')
        end
    end

    scatter(num,y_center-20,50,'r','p','filled')
    text(num*ones(1,sample),(100:-20:120-20*sample),interval_str,'FontSize',7,'HorizontalAlignment','center')
end

xlim([0.5,n_titer+0.5])
ylim([0,120]);
set(gca,'YTick',[])
set(gca,'XTick',(1:n_titer))
xlabel('vaccine combo')

%--------legend_pse--------
ax_legend=axes;
set(ax_legend,'Position',[w_b+w_left+w_jj_leg_pse,1-l_b_up-l_hum-l_hum_Rsq-l_leg_pse,w_leg,l_leg_pse])
hold on
scatter(5,8,'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_exp,'MarkerEdgeColor',c_exp,'LineWidth',1)
scatter(5,5,'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_est,'MarkerEdgeColor',c_est,'LineWidth',1)
axis off
text(6.5,8,'Experiment','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
text(6.5,5,'Model','HorizontalAlignment','left','VerticalAlignment','middle','FontSize',8)
xlim([4,10])
ylim([4,11])
%--------legend_vaccine--------
ax_legend=axes;
set(ax_legend,'Position',[w_b+w_left+w_jj_leg,l_b_leg,w_leg,l_comb])


vaccine_inform{1}='CoronaVac';
vaccine_inform{2}='RQ3013';
vaccine_inform{3}='RQ3014';
vaccine_inform{4}='RQ3021';
vaccine_inform{5}='RQ3019';
vaccine_inform{6}='RQ3033';

vaccine_inform{7}='RQ3025';
vaccine_inform{8}='RQ3027';
vaccine_inform{9}='Measurement';
ybias=0;

y_space=0.67;
y_white=1.5;
y=linspace(30,5,9);
hold on
scatter(10*ones(1,1),y(1),30,varient_color(1,:),'filled')
scatter(10*ones(1,5),y(2:6),80,varient_color(2:6,:),'s','filled')

scatter(10*ones(1,2),y(7:8)+y_space,80,varient_color(5:6,:),'s','filled')
scatter(10*ones(1,2),y(7:8)-y_space,80,varient_color([2,2],:),'s','filled')


scatter(10*ones(1,2),y(7:8)+y_white,80,'s','filled','MarkerFaceColor','w','MarkerEdgeColor','none')
scatter(10*ones(1,2),y(7:8)-y_white,80,'s','filled','MarkerFaceColor','w','MarkerEdgeColor','none')
scatter(10,y(9),50,'r','p','filled')


text(10*ones(1,9)+3,linspace(30,5,9)+ybias,vaccine_inform(1:9),'FontWeight','normal','FontSize',8)


%axis equal
xlim([10,20])
ylim([3,31])

axis off

%------Rsquare--------%
ax_Rsq=axes;
set(ax_Rsq,'Position',[1-w_b_r-w_Rsq,l_b+l_Rsq+l_jj_Rsq,w_Rsq,l_Rsq])

draw_R_square(ax_Rsq,batch_titer_lgFC,batch_titer_est_lgFC)

%------error distribution-----
ax_disb=axes;
set(ax_disb,'Position',[1-w_b_r-w_Rsq,l_b,w_Rsq,l_Rsq])

his_color=[51,186,237]/255;

lgFC_error=abs(batch_titer_est_lgFC-batch_titer_lgFC);
sum_lgFC_error=sum(lgFC_error(:),'omitnan');
mean_lgFC_error=mean(lgFC_error(:),'omitnan');
median_lgFC_error=median(lgFC_error(:),'omitnan');
hold on
histogram(lgFC_error(:),'Normalization','probability','FaceColor',his_color,'NumBins',4)
pd=fitdist(lgFC_error(:),'hn');
x_values=(0:0.05:1.2*max(lgFC_error(:)));
%pd.sigma=0.68;
y=pdf(pd,x_values);
plot(x_values,0.45*y,'r','LineWidth',1.5)

ylabel('Frequency')
xlabel('Absolute error of log_2(titer FC)')


%---------text------------
ax_tex=axes;
letter_size=10;
set(ax_tex,'Position',[0,0,1,1])
axis off
text(w_b/2,1-l_b_up/2,'a','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')

text(w_b/2,1-l_b_up-l_hum-l_hum_Rsq/2,'b','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
text(w_b+w_left+w_jj_leg+w_leg+w_jj_leg_Rsq/2,1-l_b_up-l_hum-l_hum_Rsq/2,'c','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')
text(w_b+w_left+w_jj_leg+w_leg+w_jj_leg_Rsq/2,1-l_b_up-l_hum-l_hum_Rsq-l_Rsq-l_jj_Rsq/2,'d','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',letter_size,'FontWeight','normal')