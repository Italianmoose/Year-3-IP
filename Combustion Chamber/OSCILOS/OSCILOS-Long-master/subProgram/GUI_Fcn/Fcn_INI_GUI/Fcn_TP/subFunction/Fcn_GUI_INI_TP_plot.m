function Fcn_GUI_INI_TP_plot(hAxes,handles)
% This function is used to plot the thermal propreties and mean flow along
% the combustor
% currently, the properties following are avaliable:
% 1. Mean flow veloicty
% 2. Mean temperature
%
% first created: 2014-12-03
% last modified: 2015-06-03
%
global CI
%
pop_plot = get(handles.pop_plot,'Value');           % get the value of popup menu
%------------------------------------
cla(hAxes)                                          % clear the axes
axes(hAxes)
hold on
%
N = length(CI.CD.SectionIndex);
x_plots(1,1:N-1) = CI.CD.x_sample(1:N-1);
x_plots(2,1:N-1) = CI.CD.x_sample(2:N);
%
switch pop_plot
    case 1  % mean velocity
        for ss = 1:N-1
            switch CI.CD.TubeIndex(ss)
                case 0
                    y_plots(1:2,ss) = CI.TP.u_mean(1,ss);
                case {1,2}
                    y_plots(1:2,ss) = NaN;
            end   
        end
        ylabel(hAxes,'$\overline{u}$ [m/s]','Color','k','Interpreter','LaTex','FontSize',handles.FontSize(1));
    case 2 % mean temperature
        for ss = 1:N-1
            switch CI.CD.TubeIndex(ss)
                case 0
                    y_plots(1:2,ss) = CI.TP.T_mean(1,ss);
                case {1,2}
                    y_plots(1:2,ss) = NaN;
            end  
        end
        ylabel(hAxes,'$\overline{T}$ [K]','Color','k','Interpreter','LaTex','FontSize',handles.FontSize(1));
end
%
for ss = 1:N-1
    ColorUDF{ss} = 'b';             % color of the line
end
if CI.CD.isHA == 1
    for ss = CI.CD.indexHA(1):N-1
        ColorUDF{ss} = 'r';         % after the first heat addition interface, the color of plotted lines are set to red
    end
end
% if CI.CD.isLiner == 1
%    switch pop_plot
%         case 1  
%             y_plots(2,indexLiner) = CI.TP.u_mean(1,indexLiner+1);
%         case 2
%             y_plots(2,indexLiner) = CI.TP.T_mean(1,indexLiner+1);
%    end
% end

for ss = 1:N-1
    plot(hAxes,[x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],...
        'color',ColorUDF{ss},'linewidth',2,'linestyle','-');
end
if CI.CD.isHA == 1
    for ss = 1:length(CI.CD.indexHA)
        plot(hAxes, [x_plots(1,CI.CD.indexHA(ss)),    x_plots(1,CI.CD.indexHA(ss))],...
                    [y_plots(1,CI.CD.indexHA(ss)-1),  y_plots(1,CI.CD.indexHA(ss))],...
                    'color','g','linewidth',2,'linestyle','--');
    end
end
% ------------
switch pop_plot
    case 1  % mean velocity
        for ss = 1:N-1
            a = ss
            switch CI.CD.TubeIndex(ss)
                case 0
                    % Nothing happens
                case {1,2}
                    if ss == N-1
                        y_plots(1:2,ss) = CI.TP.u_mean(1,ss);
                    else
                        y_plots(1:2,ss) = CI.TP.u_mean(1,ss:ss+1);
                    end
                    plot(hAxes,[x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],...
                        'color',ColorUDF{ss},'linewidth',2,'linestyle','-');

            end   
        end
    case 2 % mean temperature
        for ss = 1:N-1
            a = ss
            switch CI.CD.TubeIndex(ss)
                case 0
                    % Nothing happens
                case {1,2}
                    if ss == N-1
                        y_plots(1:2,ss) = CI.TP.T_mean(1,ss);
                    else
                        y_plots(1:2,ss) = CI.TP.T_mean(1,ss:ss+1);
                    end
                    plot(hAxes,[x_plots(1,ss),x_plots(2,ss)],[y_plots(1,ss),y_plots(2,ss)],...
                        'color',ColorUDF{ss},'linewidth',2,'linestyle','-');

            end   
        end
        
end
        

% ------------------------------------------------------------------------
%
set(hAxes,'YColor','k','Box','on','ygrid','on','xgrid','on');
set(hAxes,'FontName','Helvetica','FontSize',handles.FontSize(1),'LineWidth',1)
xlabel(hAxes,'$x$ [m]','Color','k','Interpreter','LaTex','FontSize',handles.FontSize(1));
set(hAxes,  'xlim',[CI.CD.x_sample(1), CI.CD.x_sample(end)],...
            'YAxisLocation','left','Color','w');
yvalue_max  = max(max(y_plots));
yvalue_min  = min(min(y_plots));  
ymax        = yvalue_max+round((yvalue_max-yvalue_min).*10)./50+eps;
ymin        = yvalue_min-round((yvalue_max-yvalue_min).*10)./50-eps;
if ymax<=ymin
    ymax = ymax+0.1*mean(mean(y_plots));
    ymin = ymin-0.1*mean(mean(y_plots));
end
set(hAxes,'ylim',[ymin ymax])
%
% -------------------------------------------------------------------------
hold off    

% --------------------------------end--------------------------------------