function pl_area_cut(x,y,Xtrunc,colors,edge_alpha,axisYLim,XLab,YLab,XTick_width,Lable_off)

    ax=gca;
    ax_width=1;
    % 找出截断区间的索引
    cut_indices = (x >= Xtrunc(1)) & x < (Xtrunc(2));

    % 创建截断后的数据
    x1 = x(~cut_indices & (x < Xtrunc(1)));
    y1 = y(:,~cut_indices & (x < Xtrunc(1)));
    x2 = x(~cut_indices & (x > Xtrunc(2)));
    y2 = y(:,~cut_indices & (x > Xtrunc(2)));

    % 绘制area图

    a=area(ax,x1,y1','EdgeAlpha',edge_alpha);
    hold(ax,'on');
    plot(ax,[min(x),Xtrunc(1)],[axisYLim(2),axisYLim(2)],'k','LineWidth',ax_width)
    ax.LineWidth=ax_width;
    for k=(1:length(a))
        a(k).FaceColor=colors(k,:);
    end
    set(ax, 'TickLength', [0.02, 0.025]);


    % box(ax,'off')
    ax.XAxisLocation='bottom';
    ax.YAxisLocation='left';

    axisPos=ax.Position;
    axisXLim=[min(x),max(x)];

    axisXScale=diff(axisXLim);
    axisYScale=diff(axisYLim);

    truncRatio=1/20; %trunc占比

    ax2=copyAxes(ax);
    
    ax2.Position=ax.Position;
    
    ax2.YTickLabels=[];
    ax2.YColor='none';


    ax.Position(3)=axisPos(3)*(1-truncRatio)/(axisXScale-diff(Xtrunc))*(Xtrunc(1)-axisXLim(1));
    ax2.Position(1)=axisPos(1)+ax.Position(3)+axisPos(3)*truncRatio;
    ax2.Position(3)=axisPos(3)*(1-truncRatio)/(axisXScale-diff(Xtrunc))*(axisXLim(2)-Xtrunc(2));

%     % 链接轴范围变动
%     linkaxes([ax,ax2],'y')
[1,1].*(ax.Position(1)+ax.Position(3))
[ax2.Position(2),ax2.Position(2)+ax2.Position(4)]
    annotation('line',[1,1].*(ax.Position(1)+ax.Position(3)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],'LineStyle',':','LineWidth',ax.LineWidth,'Color',ax.XColor);
    annotation('line',[1,1].*(ax2.Position(1)),[ax2.Position(2),ax2.Position(2)+ax2.Position(4)],'LineStyle',':','LineWidth',ax.LineWidth,'Color',ax.XColor);

    createSlash([ax.Position(1)+ax.Position(3)-.2,ax.Position(2)-.2,.4,.4])
    createSlash([ax2.Position(1)-.2,ax.Position(2)-.2,.4,.4])
    createSlash([ax.Position(1)+ax.Position(3)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4])
    createSlash([ax2.Position(1)-.2,ax.Position(2)+ax.Position(4)-.2,.4,.4])

    ax.XLim=[min(x),Xtrunc(1)];
    ax2.XLim=[Xtrunc(2),max(x)];
    ax.YLim=axisYLim;
    ax2.YLim=axisYLim;

    ax2.Position(4)=ax.Position(4);
    ax2.Position(2)=ax.Position(2);
    
    
    b=area(ax2,x2,y2','EdgeAlpha',edge_alpha);
    for k=(1:length(b))
        b(k).FaceColor=colors(k,:);
    end
    plot(ax2,[Xtrunc(1),max(x)],[axisYLim(2),axisYLim(2)],'k','LineWidth',ax_width)
    plot(ax2,[max(x),max(x)],[axisYLim(1),axisYLim(2)],'k','LineWidth',ax_width)
    ax2.LineWidth=ax_width;
    set(ax2, 'TickLength', [0.02, 0.025]);
    
    box(ax,'off')
    box(ax2,'off')
    ylabel(ax,YLab)
    
    ax.XTick=(min(x):XTick_width:Xtrunc(1));
    ax2.XTick=(Xtrunc(2):XTick_width:max(x));
    
    ax.FontSize=6;
    ax2.FontSize=6;
    ax.YLabel.FontSize=8;
    ax2.YLabel.FontSize=8;
    set(ax,'TickDir','out')
    set(ax2,'TickDir','out')
    if Lable_off
        ax.XTickLabel=[];
        ax2.XTickLabel=[];

    end
    
    
    function newAX=copyAxes(ax)
        axStruct=get(ax);
        fNames=fieldnames(axStruct);
        newAX=axes('Parent',ax.Parent);

        coeList={'CurrentPoint','XAxis','YAxis','ZAxis','BeingDeleted',...
            'TightInset','NextSeriesIndex','Children','Type','Legend'};
        for n=1:length(coeList)
            coePos=strcmp(fNames,coeList{n});
            fNames(coePos)=[];
        end

        for n=1:length(fNames)
            newAX.(fNames{n})=ax.(fNames{n});
        end

        copyobj(ax.Children,newAX)
    end

    function createSlash(pos)
        anno=annotation('textbox');
        anno.String='/';
        anno.LineStyle='none';
        anno.FontSize=6;
        anno.Position=pos;
        anno.FitBoxToText='on';
        anno.VerticalAlignment='middle';
        anno.HorizontalAlignment='center';
    end
end
