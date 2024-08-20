clc,clear,close all
%load('IPD_pre_amount_distance_mRNA_boost.mat')
load('fig6a_data.mat')
Ab_1_sample=zeros(length(pre_amount),length(distance));
Ab_2_sample=zeros(length(pre_amount),length(distance));
M_1_sample=zeros(length(pre_amount),length(distance));
M_2_sample=zeros(length(pre_amount),length(distance));
titer_2_sample=zeros(length(pre_amount),length(distance));
titer_1_sample=zeros(length(pre_amount),length(distance));
X=eye(2,2);
self_neu=ones(2,1);
vaccine_time=[10,40,400,430];
gap_sample=14;

for i=(1:length(pre_amount))
    i
    for j=(1:length(distance))
        cross_neu=2^(-distance(j));
        X(1,2)=cross_neu;
        X(2,1)=cross_neu;  
        t=t_amount_dis{i,j};
        Agt=Ag_amount_dis{i,j};
        Abt=Ab_amount_dis{i,j};
        Ft=F_amount_dis{i,j};
        M_offt=M_off_amount_dis{i,j};
        M_ont=M_on_amount_dis{i,j};
        Mt=M_offt+M_ont;
        titert_2=X(2,:)'.*Abt.*self_neu;
        titert_1=X(1,:)'.*Abt.*self_neu;
        
        t_sample=vaccine_time(end)+gap_sample;
        [Min,ind]=min(abs(t-t_sample'),[],2);
        titer_2_sample(i,j)=sum(titert_2(:,ind));
        titer_1_sample(i,j)=sum(titert_1(:,ind));
        
        Ab_1_sample(i,j)=Abt(1,ind);
        Ab_2_sample(i,j)=Abt(2,ind);
        M_1_sample(i,j)=Mt(1,ind);
        M_2_sample(i,j)=Mt(2,ind);        
    end 
end
lg2_titer_2=log2(titer_2_sample./titer_2_sample(:,1));
lg2_titer_1=log2(titer_1_sample./titer_1_sample(:,1));

[M_titer,I_titer]=min(lg2_titer_2,[],2);

num=0;
l_bound_titer_all=zeros(length(pre_amount),1);
r_bound_titer_all=zeros(length(pre_amount),1);
for check_amount=pre_amount
    num=num+1;
    I_titer_check=I_titer(pre_amount==check_amount);
    [~,l_bound_titer]=min(abs(lg2_titer_2(pre_amount==check_amount,1:I_titer_check)+1),[],2);
    [~,r_bound_titer]=min(abs(lg2_titer_2(pre_amount==check_amount,I_titer_check:end)+1),[],2);
    r_bound_titer=I_titer_check+r_bound_titer-1;
    l_bound_titer_all(num)=l_bound_titer;
    r_bound_titer_all(num)=r_bound_titer;
    
end
%%
ad_range=[0,8];
ad_rangel=-2.^(-ad_range);
distancel=-2.^(-distance);
c_50pro=[ 150, 150, 150 ]/255;

close all
figure(1)
set(gcf,'Position',[600,200,900,750])

ax_inten=axes;
set(ax_inten,'Position',[0.22,0.1,0.72,0.75])

[XX,YY]=meshgrid(distance,pre_amount);
XXl=-(2.^(-XX));
s=surface(XXl,YY,lg2_titer_2,'EdgeColor','k','FaceColor','interp','EdgeAlpha',0.2);
%s=surface(XX,YY,lg2_titer_2,'EdgeColor','k','FaceColor','interp');
view(3)
set(ax_inten,'XScale','log')
set(ax_inten,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','0.5','0.75','0.9','0.95','0.975','0.99','0.995','0.998'})
set(ax_inten,'YDir','reverse')
colormap(ax_inten,flip(othercolor('RdBu4')));
% s.CData=map;

xlabel('Antigenic distance')
xlim(ad_rangel)
ylim([1,40])
ylabel('Primary dosage \it d \rm(\mug mRNA vaccine)')
zlabel('Added neutralization against variant')

hold on
[M_titer,I_titer]=min(lg2_titer_2,[],2);



hold on
surface(ax_inten,XXl,YY,-1*ones(size(XX)),'EdgeColor','none','FaceColor',c_50pro,'FaceAlpha',0.25)


%text(ax_inten,-2^-4.4,36,'minimum protection','Color','r','VerticalAlignment','bottom','HorizontalAlignment','center')


c_amount=[0.4940 0.1840 0.5560];
c_amount_weak=[255 20 147]/255;
plot3(distancel,30*ones(size(distancel)),lg2_titer_2(pre_amount==30,:),'-','LineWidth',2,'Color',c_amount)

text(ax_inten,-1,30,'Strong primary immunization (\itd=30)','VerticalAlignment','bottom','HorizontalAlignment','right','Color',c_amount,'FontWeight','bold')

plot3(distancel,3*ones(size(distancel)),lg2_titer_2(pre_amount==3,:),'-','LineWidth',2,'Color',c_amount_weak)

text(ax_inten,-1,3,'Weak primary immunization (\itd=3)','VerticalAlignment','bottom','HorizontalAlignment','right','Color',c_amount_weak,'FontWeight','bold')
set(ax_inten,'TickDir','out')

z_bound=-4.5;
set(ax_inten,'ZTick',[-3,-2,-1,0],'ZTickLabel',{'12.5%','25%','50%','100%'},'ZLim',[z_bound,0])



x_q=(0:0.1:10);
y_q=(1:1:50);
[XX,YY]=meshgrid(x_q,y_q);
XXl=-(2.^(-XX));

%-----plot left bound------
xx_left=distance(l_bound_titer_all);
yy_left=pre_amount;
xx_left_q = interp1(yy_left,xx_left,y_q);

xx_left_ql=-2.^(-xx_left_q);


xx_left_all=distance(l_bound_titer_all);

xx_left=distance(l_bound_titer_all(1:2:40));
yy = smooth(xx_left,pre_amount(1:2:40),'sgolay');

xx_left=[xx_left,xx_left(end)];
xx_leftl=-2.^(-xx_left);
yy=[yy;40];

plot3(xx_leftl,yy,z_bound*ones(size(yy)),'--','LineWidth',2,'Color',c_50pro)

z_50=log2(0.5);

%plot3(xx_leftl,yy,z_50*ones(size(yy)),'--','LineWidth',2,'Color',c_50pro)

%-----plot right bound------


xx_right=distance(r_bound_titer_all);
yy_right=pre_amount;
xx_right_q = interp1(yy_right,xx_right,y_q);

xx_right_ql=-2.^(-xx_right_q);

xx_right=distance(r_bound_titer_all);
xx_rightl=-2.^(-xx_right);
yy = smooth(xx_rightl,pre_amount);

plot3(xx_rightl,yy,z_bound*ones(size(yy)),'--','LineWidth',2,'Color',c_50pro)
%plot3(xx_right_ql,y_q,z_bound*ones(size(y_q)),'--','LineWidth',2,'Color',c_50pro)
%plot3(xx_rightl,yy,z_50*ones(size(yy)),'--','LineWidth',2,'Color',c_50pro)





ind_X=repmat((1:101),50,1);
ind_pro=ind_X<l_bound_titer_all;
ind_brk=ind_X>r_bound_titer_all;
ind_pit=ind_X>=l_bound_titer_all|ind_X<=r_bound_titer_all;

min_pro=min(lg2_titer_2(ind_pro));
max_pro=max(lg2_titer_2(ind_pro));

min_pit=min(lg2_titer_2(ind_pit));
max_pit=max(lg2_titer_2(ind_pit));

min_brk=min(lg2_titer_2(ind_brk));
max_brk=max(lg2_titer_2(ind_brk));


ind_X=repmat(x_q,length(y_q),1);
ind_pro=ind_X<xx_left_q';
ind_brk=ind_X>xx_right_q';
ind_pit=ind_X>=xx_left_q'|ind_X<=xx_right_q';


CC=0.75*ones([size(XX),3]);

c_pro=[linspace(230,139,256)',linspace(225,189,256)',linspace(226,164,256)']/256;
c_pit=[linspace(255,216,256)',linspace(230,121,256)',linspace(218,127,256)']/256;
c_brk=[linspace(230,118,256)',linspace(256,157,256)',linspace(256,234,256)']/256;


for i=(1:length(y_q))
    for j=(1:length(x_q))
        val=lg2_titer_2(ceil(i/1),ceil(j/1));
        if ind_pro(i,j)
            norm_value=(val-min_pro)/(max_pro-min_pro);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_pro(c_ind,:);
            
        elseif ind_brk(i,j)
            norm_value=(val-min_brk)/(max_brk-min_brk);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_brk(c_ind,:);
        else
            norm_value=1-(val-min_pit)/(max_pit-min_pit);
            c_ind=round(256*norm_value);
            if c_ind<=0
                c_ind=1;
            end
            CC(i,j,:)=c_pit(c_ind,:);
        end
    end
end

surface(ax_inten,XXl,YY,z_bound*ones(size(XX)),CC,'EdgeColor','none','FaceLighting','none','FaceAlpha',0.45)


y_tex_inten=27.5;
z_tex=-3.75;
fsize=12;
x_adj=-0.8;

text(ax_inten,-2^(-distance(l_bound_titer)/2+x_adj),y_tex_inten,z_tex,{'Immune-Imprinting- ';'Protection Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.5,0],'FontSize',fsize,'FontWeight','bold')


text(ax_inten,-2^(-(distance(l_bound_titer)+distance(r_bound_titer))/2+x_adj),y_tex_inten,z_tex,{'Immune-Imprinting-';'Pitfall Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.7,0,0],'FontSize',fsize,'FontWeight','bold')


text(ax_inten,-2^(-(ad_range(2)+distance(r_bound_titer))/2+x_adj),y_tex_inten,z_tex,{'Immune-Imprinting-';'Breakthrough Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.5],'FontSize',fsize,'FontWeight','bold')

%text(ax_inten,-0.0025,40,-0.8,'50% escape','Color',0.8*c_50pro,'VerticalAlignment','bottom','HorizontalAlignment','center')
view(ax_inten,[25.1188021427944 21.9067415730338]);