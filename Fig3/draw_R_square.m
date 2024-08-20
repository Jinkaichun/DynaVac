function draw_R_square(ax,batch_titer_lgFC,batch_titer_est_lgFC,batch_dataset)
alpha=0.5;
hold on

c_yd=[0,0,255]/255;
c_cyl=[0,150,0]/255;

scatter(ax,batch_titer_lgFC(batch_dataset==1),batch_titer_est_lgFC(batch_dataset==1),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_yd)
scatter(ax,batch_titer_lgFC(batch_dataset==0),batch_titer_est_lgFC(batch_dataset==0),'filled','MarkerFaceAlpha',alpha,'MarkerFaceColor',c_cyl)



min_FC=min([batch_titer_lgFC(:),batch_titer_est_lgFC(:)],[],2,'includenan');
max_FC=max([batch_titer_lgFC(:),batch_titer_est_lgFC(:)],[],2,'includenan');

min_FC=min(min_FC);
max_FC=max(max_FC);


lgFC_width=[min_FC-0.5,max_FC+0.5];

plot(ax,lgFC_width,lgFC_width,'--r','LineWidth',2)
set(ax,'XTick',(floor(min_FC):ceil(max_FC)),'YTick',(floor(min_FC):ceil(max_FC)),'TickDir','out')

mean_lgFC=mean(batch_titer_lgFC(:),'omitnan');
TSS=sum((batch_titer_lgFC(:)-mean_lgFC).^2,'omitnan');
ESS=sum((batch_titer_lgFC(:)-batch_titer_est_lgFC(:)).^2,'omitnan');
R=1-ESS/TSS;
n_data=length(batch_titer_est_lgFC(:))-sum(isnan(batch_titer_lgFC(:)));

fsize=8;
text(ax,0.05,0.82,['R^2=',num2str(round(R,3))],'Units','normalized','HorizontalAlignment','left','FontSize',fsize)
text(ax,0.05,0.9,['n=',num2str(n_data)],'Units','normalized','HorizontalAlignment','left','FontSize',fsize)
xlabel(ax,'log_2 FC(neutralization) (Experiment) ')
ylabel(ax,'log_2 FC(neutralization) (Model) ')
xlim(ax,lgFC_width)
ylim(ax,lgFC_width)
end