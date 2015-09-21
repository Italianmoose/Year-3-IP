function Fcn_CD_plot_with_dampers(hAxes,handles,indexLegend)
% This function is used to plot the schematic view of the combustor
% The input are the tag of axes and handles of current gui window
% first created: 2014-12-03
%
global CI
hFontsize2  = handles.FontSize(2);
poshAxes    = get(handles.axes1,'position');
WallLineWidth = 2;
% -------------------------------------------------------------------------
%
x_sample =      CI.CD.x_sample;                
r_sample =      CI.CD.r_sample;
SectionIndex=   CI.CD.SectionIndex;

% Check if damper exists. Initialize the damper setting if damper does not
% exist. This is not needed in GUI_INI_PD, but should be included in
% GUI_INI_CD.
idExistHR=any(strcmp('NUM_HR',fieldnames(CI.CD)));
idExistLiner=any(strcmp('NUM_Liner',fieldnames(CI.CD)));
if idExistHR==0
CI.CD.NUM_HR      = 0;
CI.CD.HR_config   = [];
else
CI.CD.NUM_HR      = CI.CD.NUM_HR;
CI.CD.HR_config   = CI.CD.HR_config;
end

if idExistLiner==0
CI.CD.NUM_Liner      = 0;
CI.CD.Liner_config   = [];
else
CI.CD.NUM_Liner      = CI.CD.NUM_Liner;
CI.CD.Liner_config   = CI.CD.Liner_config;
end
% Define Liner_config for Liner plotting
NUM_Liner=CI.CD.NUM_Liner;
if NUM_Liner==0
    Liner_config=[];
else
    Liner_config=CI.CD.Liner_config;
end


%-------------------------------------
W           = abs(x_sample(end) - x_sample(1));             % Length of the combustor
H           = 3.5*max(r_sample);                              % Diameter of the combustor
plot_ratio  = 1.1;                                          % Ratio of the axes limit to the combustor dimension
axes_W      = plot_ratio*W;                                 % axes x width
axes_H      = 2*plot_ratio*H;                               % axes y width        
x_min       = x_sample(1) - (axes_W-W)./2;                  % axes x min
y_min       = -axes_H./2;                                   % axes y min
%--------------------------------------
cla(hAxes);                             % clear the current axes 
axis(hAxes);
hold on
% -------------------------------------
% These are used for the legend
xMax = max(get(hAxes,'xlim'));
yMax = max(get(hAxes,'ylim'));
plot(hAxes,xMax,yMax,'-b','linewidth',3)
plot(hAxes,xMax,yMax,'-g','linewidth',3)
plot(hAxes,xMax,yMax,'-r','linewidth',3)
plot(hAxes,xMax,yMax,'-m','linewidth',3)
%--------------------------------------
% plot the approximate profile of the combustor which consisting of several
% sections 



NUM_Liner_note=1;
NUM_HR_note=1;
% Configure the plotting scale of the HRs
HR_neck_width=axes_W/150;
HR_neck_length=0.15*max(r_sample);
HR_cavity_width=axes_W/50;
HR_cavity_length=0.8*max(r_sample);
for s=1:length(x_sample)-1
    % Plot the outline of the combustor
        x_plot1=1000*[x_sample(s),x_sample(s+1)];
        x_plot2=1000*[x_sample(s),x_sample(s+1)-HR_neck_width/2];
        x_plot3=1000*[x_sample(s)+HR_neck_width/2,x_sample(s+1)-HR_neck_width/2];
        x_plot4=1000*[x_sample(s)+HR_neck_width/2,x_sample(s+1)];
        x_plot5=1000*[x_sample(s+1),x_sample(s+1)];
        
        y_plot1=1000*[r_sample(s),r_sample(s)];
        y_plot2=1000*[-r_sample(s),-r_sample(s)];
        y_plot3=1000*[r_sample(s),r_sample(s+1)];
        y_plot4=1000*[-r_sample(s),-r_sample(s+1)];
    if SectionIndex(s)==2
        if SectionIndex(s+1)==2
          plot(hAxes,x_plot3, y_plot1,'-k','linewidth',WallLineWidth);
          plot(hAxes,x_plot1, y_plot2,'-k','linewidth',WallLineWidth);
        else
                plot(hAxes,x_plot4, y_plot1,'-k','linewidth',WallLineWidth);
                plot(hAxes,x_plot1, y_plot2,'-k','linewidth',WallLineWidth);
        end
    else if SectionIndex(s)==30
                plot(hAxes,x_plot1, y_plot1,'--k','linewidth',WallLineWidth);
                plot(hAxes,x_plot1, y_plot2,'--k','linewidth',WallLineWidth);
        
    else if SectionIndex(s)==31
           if SectionIndex(s+1)==2
               plot(hAxes,x_plot2, y_plot1,'-k','linewidth',WallLineWidth);
               plot(hAxes,x_plot1, y_plot2,'-k','linewidth',WallLineWidth);
           else 
                plot(hAxes,x_plot1, y_plot1,'-k','linewidth',WallLineWidth);
                plot(hAxes,x_plot1, y_plot2,'-k','linewidth',WallLineWidth);
           end
    else if SectionIndex(s+1)==2 
          plot(hAxes,x_plot2, y_plot1,'-k','linewidth',WallLineWidth);
          plot(hAxes,x_plot1, y_plot2,'-k','linewidth',WallLineWidth);
    else
                plot(hAxes,x_plot1, y_plot1,'-k','linewidth',WallLineWidth);
                plot(hAxes,x_plot1, y_plot2,'-k','linewidth',WallLineWidth);
        end
        end
        end
    end
    
    if r_sample(s) ~=r_sample(s+1)
        plot(hAxes,x_plot5, y_plot3,'-k','linewidth',WallLineWidth);
        plot(hAxes,x_plot5, y_plot4,'-k','linewidth',WallLineWidth);
    else
    end
    
    % plot the HR
    if SectionIndex(s)==2
        x_plot=1000*[x_sample(s)-HR_neck_width/2, x_sample(s)-HR_neck_width/2, x_sample(s)-HR_cavity_width/2,  x_sample(s)-HR_cavity_width/2,               x_sample(s)+HR_cavity_width/2,               x_sample(s)+HR_cavity_width/2, x_sample(s)+HR_neck_width/2, x_sample(s)+HR_neck_width/2];
        y_plot=1000*[r_sample(s),                 r_sample(s)+HR_neck_length,  r_sample(s)+HR_neck_length,    r_sample(s)+HR_neck_length+HR_cavity_length,  r_sample(s)+HR_neck_length+HR_cavity_length, r_sample(s)+HR_neck_length,    r_sample(s)+HR_neck_length,  r_sample(s)];
        plot(hAxes,x_plot, y_plot,'-k','linewidth',WallLineWidth);
        text(1000*x_sample(s), -1000*1.5*max(r_sample),...
             ['HR',num2str(NUM_HR_note)],'FontSize',9,...
             'interpreter','latex','HorizontalAlignment','center');
        NUM_HR_note=NUM_HR_note+1; 
    end
    % plot the Liner
        if SectionIndex(s)==30
            if Liner_config(NUM_Liner_note,4)==0
                    if Liner_config(NUM_Liner_note,5)==0
                            x_plot=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plot=1000*[r_sample(s),    r_sample(s)*3,  r_sample(s)*3, r_sample(s)];
                            plot(hAxes,x_plot, y_plot,'-k','linewidth',WallLineWidth);
                            x_plotn=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plotn=1000*[-r_sample(s),-r_sample(s)*3,-r_sample(s)*3,-r_sample(s)];
                            plot(hAxes,x_plotn, y_plotn,'-k','linewidth',WallLineWidth);
                            text(1000*(x_sample(s)+x_sample(s+1))/2, -1000*3.5*max(r_sample),...
                                ['Liner',num2str(NUM_Liner_note)],'FontSize',9,...
                                'HorizontalAlignment','center');
                            NUM_Liner_note=NUM_Liner_note+1;
                    else
                            x_plot=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plot=1000*[r_sample(s),    r_sample(s)*1.2,  r_sample(s)*1.2, r_sample(s)];
                            plot(hAxes,x_plot, y_plot,'-k','linewidth',WallLineWidth);
                            x_plotn=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plotn=1000*[-r_sample(s),-r_sample(s)*1.2,-r_sample(s)*1.2, -r_sample(s)];
                            plot(hAxes,x_plotn, y_plotn,'-k','linewidth',WallLineWidth);
                            text(1000*(x_sample(s)+x_sample(s+1))/2, -1000*3.5*max(r_sample),...
                                ['Liner',num2str(NUM_Liner_note)],'FontSize',9,...
                                'HorizontalAlignment','center');
                            NUM_Liner_note=NUM_Liner_note+1;
                    end     
            else
                    if Liner_config(NUM_Liner_note,5)==0
                            x_plot=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plot=1000*[r_sample(s),    r_sample(s)*3,  r_sample(s)*3, r_sample(s)];
                            x_plot2=1000*[x_sample(s), x_sample(s+1)];
                            y_plot2=1000*[r_sample(s)*1.3,  r_sample(s)*1.3];                            
                            plot(hAxes,x_plot, y_plot,'-k','linewidth',WallLineWidth);
                            plot(hAxes,x_plot2, y_plot2,'--k','linewidth',WallLineWidth);
                            x_plotn=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plotn=1000*[-r_sample(s),    -r_sample(s)*3,  -r_sample(s)*3, -r_sample(s)];
                            x_plot2n=1000*[x_sample(s), x_sample(s+1)];
                            y_plot2n=1000*[-r_sample(s)*1.3,  -r_sample(s)*1.3];                            
                            plot(hAxes,x_plotn, y_plotn,'-k','linewidth',WallLineWidth);
                            plot(hAxes,x_plot2n, y_plot2n,'--k','linewidth',WallLineWidth);
                            text(1000*(x_sample(s)+x_sample(s+1))/2, -1000*3.5*max(r_sample),...
                                ['Liner',num2str(NUM_Liner_note)],'FontSize',9,...
                                'HorizontalAlignment','center');
                            NUM_Liner_note=NUM_Liner_note+1;
                    else
                            x_plot=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plot=1000*[r_sample(s),    r_sample(s)*1.5,  r_sample(s)*1.5, r_sample(s)];
                            x_plot2=1000*[x_sample(s), x_sample(s+1)];
                            y_plot2=1000*[r_sample(s)*1.3,  r_sample(s)*1.3];                            
                            plot(hAxes,x_plot, y_plot,'-k','linewidth',WallLineWidth);
                            plot(hAxes,x_plot2, y_plot2,'--k','linewidth',WallLineWidth);
                            x_plotn=1000*[x_sample(s), x_sample(s), x_sample(s+1),  x_sample(s+1)];
                            y_plotn=1000*[-r_sample(s),    -r_sample(s)*1.5,  -r_sample(s)*1.5, -r_sample(s)];
                            x_plot2n=1000*[x_sample(s), x_sample(s+1)];
                            y_plot2n=1000*[-r_sample(s)*1.3,  -r_sample(s)*1.3];                            
                            plot(hAxes,x_plotn, y_plotn,'-k','linewidth',WallLineWidth);
                            plot(hAxes,x_plot2n, y_plot2n,'--k','linewidth',WallLineWidth);
                            text(1000*(x_sample(s)+x_sample(s+1))/2, -1000*3.5*max(r_sample),...
                                ['Liner',num2str(NUM_Liner_note)],'FontSize',9,...
                                'HorizontalAlignment','center');
                            NUM_Liner_note=NUM_Liner_note+1;
                    end
            end
        end 

end
% ---



for s = 1:length(CI.CD.SectionIndex)
    switch CI.CD.SectionIndex(s)     % define the color of different sectional interfaces
        case 0          % just interface
            indexColor = 'k';
            indexLinestyle  = '-';
            indexLineWidth  = 1;
        case 2          % HR
            indexLinestyle = 'none';
        case 10         % with heat addition but perturbation
            indexColor = 'm';
            indexLinestyle  = '-';
            indexLineWidth  = 3;
        case 11         % with heat addition and heat perturbations
            indexColor = 'r';
            indexLinestyle  = '-';
            indexLineWidth  = 3;
    end 
    if s == 1
        indexColor = 'b';   % inlet
        indexLinestyle  = '-';
        indexLineWidth  = 3;
    elseif s == length(CI.CD.SectionIndex)
        indexColor = 'g';   % outlet
        indexLinestyle  = '-';
        indexLineWidth  = 3;
    end
    if CI.CD.TubeIndex(s) == 1
        indexLinestyle  = 'none';
    end
    if CI.CD.TubeIndex(s) == 2
        indexLinestyle  = 'none';
    end
    plot(hAxes,1e3*[x_sample(s),x_sample(s)],-1e3*[-r_sample(s),r_sample(s)],...
        'linestyle',indexLinestyle,...
        'color',indexColor,'linewidth',indexLineWidth);
end
%
% in case CI.CD.TubeIndex(s) == 1
diffTubeIndex = diff(CI.CD.TubeIndex);
indexVarTubeIndex = find(diffTubeIndex~=0);
if ~isempty(indexVarTubeIndex)
    for k = 1:length(indexVarTubeIndex)
        s = indexVarTubeIndex(k)+1;
%         if CI.CD.TubeIndex(s-1)==2 || CI.CD.TubeIndex(s)==2
%             % Do nothing
%         else
        plot(hAxes,1e3*[x_sample(s),x_sample(s)],-1e3*[-r_sample(s),r_sample(s)],...
        'linestyle','-',...
        'color','k','linewidth',1);
%         end
    end
end
set(hAxes,'xlim',1000*[x_min, x_min+axes_W]);
set(hAxes,'box','on','linewidth',0.5,'gridlinestyle','-.');
set(hAxes,'xgrid','on','ygrid','on');
set(hAxes,'ylim',1000*[y_min, y_min+axes_H]);
xlabel(hAxes,'x~ [mm]','Color','k','Interpreter','LaTex');
ylabel(hAxes,'r~ [mm]','Color','k','Interpreter','LaTex');   
%
switch indexLegend
    case 1
    newline = char(10);
    legend1 = ['with mean heat', newline, 'addition and', newline, 'heat perturbations'];
    legend2 = ['with mean heat', newline, 'addition but no', newline, 'heat perturbation'];
    hlegend = legend(hAxes,...
                            'inlet','outlet',...
                            legend1,...
                            legend2);
    set(hlegend,'fontsize',hFontsize2,'location','northeastoutside');
    otherwise
end
set(handles.axes1,      'units', 'points',...
                        'Fontunits','points',...
                        'position',poshAxes); 
%
% ------------------------------end----------------------------------------              