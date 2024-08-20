clc,clear,close all
load('fig5_data.mat')
%%
clc,close all
ad_range=[0,8];
ad_rangel=-2.^(-ad_range);
distance=(0:0.1:8);
amount=(1:50);
%-----------3d-------------

figure(2)
set(gcf,'Position',[400,200,850,750])
ax_inten=axes;
set(ax_inten,'Position',[0.1,0.1,0.75,0.75])

[XX,YY]=meshgrid(distance,amount);
XXl=-(2.^(-XX));
s=surface(XXl,YY,titer_2_sample,'EdgeColor','k','FaceColor','interp','EdgeAlpha',0.2);
%s=surface(XX,YY,lg2_titer_2,'EdgeColor','k','FaceColor','interp');
view(3)
set(ax_inten,'XScale','log')
set(ax_inten,'XTick',-[1,0.5,0.25,0.1,0.05,0.025,0.01,0.005,0.002],'XTickLabel',{'0','0.5','0.75','0.9','0.95','0.975','0.99','0.995','0.998'})
set(ax_inten,'YDir','normal')
colormap(ax_inten,flip(othercolor('RdBu4')));
% s.CData=map;

xlabel('Antigenic distance')
xlim(ad_rangel)
ylim([1,40])
ylabel('Booster dosage (\mug mRNA vaccine)')
zlabel('Added neutralization against variant')

z_bound=-1200;
set(ax_inten,'ZTick',(0:500:3500),'ZTickLabel',(0:500:3500),'ZLim',[z_bound,3500])



CC=0.75*ones([size(XX),3]);


titer_2_sample_check=titer_2_sample(amount==30,:);


ind_pro=distance<2.1;
ind_brk=distance>=5.5;
ind_pit=distance<5.5 &distance>=2.1;

% ind_pro=repmat(distance<2.1,50,1);
% ind_brk=repmat(distance>5.5,50,1);
% ind_pit=repmat(distance<=5.5 &distance>=2.1,50,1);


min_pro=min(titer_2_sample_check(ind_pro));
max_pro=max(titer_2_sample_check(ind_pro));

min_pit=min(titer_2_sample_check(ind_pit));
max_pit=max(titer_2_sample_check(ind_pit));

min_brk=min(titer_2_sample_check(ind_brk));
max_brk=max(titer_2_sample_check(ind_brk));

c_pro=othercolor('Greens3');
c_pit=othercolor('Reds3');
c_brk=othercolor('Bu_10');

for i=(1:length(distance))
    val=titer_2_sample_check(i);
    if ind_pro(i)
        norm_value=(val-min_pro)/(max_pro-min_pro);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC(:,i,:)=reshape(repmat(c_pro(c_ind,:),50,1),50,1,3);

    elseif ind_brk(i)
        norm_value=(val-min_brk)/(max_brk-min_brk);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC(:,i,:)=reshape(repmat(c_brk(c_ind,:),50,1),50,1,3);
    else
        norm_value=1-(val-min_pit)/(max_pit-min_pit);
        c_ind=round(256*norm_value);
        if c_ind<=0
            c_ind=1;
        end
        CC(:,i,:)=reshape(repmat(c_pit(c_ind,:),50,1),50,1,3);
    end
end


surface(ax_inten,XXl,YY,z_bound*ones(size(XX)),CC,'EdgeColor','none','FaceLighting','none','FaceAlpha',0.5)
hold on

c_50pro=[ 110, 110, 110 ]/255;

zz=(z_bound:100:3500);
[ZZ,YY2]=meshgrid(zz,amount);

y_tex_inten=12.5;
z_tex=-900;
fsize=12;


plot3(ax_inten,(-2^(-2.1))*ones(size(ZZ)),amount,titer_2_sample(:,distance==2.1),'--','LineWidth',2,'Color',c_50pro);
plot3(ax_inten,(-2^(-5.5))*ones(size(ZZ)),amount,titer_2_sample(:,distance==5.5),'--','LineWidth',2,'Color',c_50pro);
surface(ax_inten,(-2^(-2.1))*ones(size(ZZ)),YY2,ZZ,'FaceColor',1.5*c_50pro,'FaceAlpha',0.5,'EdgeColor','none');
surface(ax_inten,(-2^(-5.5))*ones(size(ZZ)),YY2,ZZ,'FaceColor',1.5*c_50pro,'FaceAlpha',0.5,'EdgeColor','none');
text(ax_inten,-2^(-distance(22)/2),y_tex_inten,z_tex,{'Immune-Imprinting- ';'Protection Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.3,0],'FontSize',fsize,'FontWeight','bold')

text(ax_inten,-2^(-(distance(22)+distance(56))/2),y_tex_inten,z_tex,{'Immune-Imprinting-';'Pitfall Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.5,0,0],'FontSize',fsize,'FontWeight','bold')

text(ax_inten,-2^(-(ad_range(2)+distance(56))/2),y_tex_inten,z_tex,{'Immune-Imprinting-';'Breakthrough Zone'},'HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.3],'FontSize',fsize,'FontWeight','bold')


plot3((-2^(-1))*ones(size(ZZ)),amount,titer_2_sample(:,distance==1),'-','LineWidth',2,'Color',[0,0.5,0])
plot3((-2^(-4.3))*ones(size(ZZ)),amount,titer_2_sample(:,distance==4.3),'-','LineWidth',2,'Color',[0.7,0,0])
plot3((-2^(-6.7))*ones(size(ZZ)),amount,titer_2_sample(:,distance==6.7),'-','LineWidth',2,'Color',[0,0,0.5])

z_tex2=3500;
% text(ax_inten,-2^(-distance(22)/2),y_tex_inten,z_tex2,'Cross-Neutralization=50%','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0.5,0],'FontSize',10,'FontWeight','bold')
% 
% text(ax_inten,-2^(-(distance(30)+distance(56))/2),y_tex_inten,z_tex2-2000,'Cross-Neutralization=5%','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0.7,0,0],'FontSize',10,'FontWeight','bold')
% 
% text(ax_inten,-2^(-(ad_range(2)+distance(60))/2),y_tex_inten,z_tex2,'Cross-Neutralization=1%','HorizontalAlignment','center','VerticalAlignment','middle','Color',[0,0,0.5],'FontSize',10,'FontWeight','bold')

view(ax_inten,[-353.546197857206 13.59558772688]);