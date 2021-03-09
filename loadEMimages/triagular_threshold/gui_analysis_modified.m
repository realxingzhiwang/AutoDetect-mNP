function particleanalysis2_1
close all;
%%Nanoparticle Sizing and Morphology Analysis Program%%
%%by Alex Powers 
%%last edit 1/26/2015


%%%%%%%%Initilization Tasks%%%%%%%%%%%
screen = get(0, 'ScreenSize'); 
width = 1100; 
height = 500-145; 
left = screen(3)/2 - width/2; 
bottom = screen(4)/2 - height/2; 


%%%%%%%%Component Construction %%%%%%%
    f1 = figure('Position', [left-100, bottom+125, 3*(300)+27.5*4, height+25], 'Visible', 'on','Color',[0.9411 0.9411 0.9411]); %f1 contains graphs
    f2 = figure('Position', [1120, 350, 238, 380], 'Visible', 'on','Color',[0.9411 0.9411 0.9411]); 
    f3 = figure('Color',[0.9411 0.9411 0.9411], 'ToolBar', 'Figure'); 
    
    set(0, 'currentfigure', f1);
    
    %create 3 axis
    aborder = 27.5; 
    newaborder = 27; 
    awidth = 300-2*aborder;
    newawidth = 330-2*aborder; 
    abottom = 27;
    aheight = height - abottom - aborder; 
    ha1 = axes('Units', 'pixels', 'Position',[2*newaborder, abottom+15, newawidth, aheight]); 
    ha2 = axes('Units', 'pixels', 'Position',[4*newaborder+newawidth, abottom+15, newawidth, aheight]);
    ha3 = axes('Units', 'pixels', 'Position',[6*newaborder+2*newawidth, abottom+15, newawidth, aheight]);

    set(0, 'currentfigure', f2);
    
    %FITTING BUTTONS
    fittingpanel = uipanel('Title','Fitting','Position',[0.038, 0.024, 0.824, 0.261],'BorderType','etchedin','Visible','on',...
             'BackgroundColor',[0.9411 0.9411 0.9411]);
    fitlabel = uicontrol('Style', 'text', 'String', 'Range to Fit', 'FontSize',8, 'Position', [90 72 88 15]);
    
    %axis 1: checkbox for fit, limits of fit
    checkbox1 = uicontrol('Style', 'checkbox', 'String', 'Diameter','Position',[21, 49, 75, 23],'Callback', {@checkbox1_Callback} );    
    fit1min = uicontrol('Style', 'edit', 'String','2', 'Position',[90, 49, 39,18], 'Callback', {@fit1min_Callback}, 'BackgroundColor', 'white'); 
    fit1max = uicontrol('Style', 'edit', 'String','15', 'Position',[140, 49 39 18], 'Callback', {@fit1max_Callback}, 'BackgroundColor', 'white');
   
    %axis 2: checkbox for fit, limits of fit
    shift = awidth + 2*aborder; 
    checkbox2 = uicontrol('Style', 'checkbox', 'String', 'Feret', 'Position',[22, 23, 75, 23],'Callback', {@checkbox2_Callback} );
    fit2min = uicontrol('Style', 'edit', 'String','0', 'Position',[90 23 38 18], 'Callback', {@fit2min_Callback}, 'BackgroundColor', 'white');
    fit2max = uicontrol('Style', 'edit', 'String','100', 'Position',[140 23 38 18], 'Callback', {@fit2max_Callback}, 'BackgroundColor', 'white');
    
    %FILTERING BUTTONS
    %roundness cuttoff 
    filterpanel = uipanel('Title','Filtering','Position',[0.034 0.3 0.836 0.482],'BorderType','etchedin','Visible','on',...
             'BackgroundColor',[0.9411 0.9411 0.9411]);
    
    filtlabel1 = uicontrol('Style', 'text', 'String', 'Min.', 'FontSize',9, 'Position', [93 260 39 18]);
    filtlabel2 = uicontrol('Style', 'text', 'String', 'Max.', 'FontSize',9, 'Position', [150 260 39 18]);
    
    rndminbutton = uicontrol('Style', 'edit', 'String','0', 'Position',...
        [91,163,42,24], 'Callback', {@rdcutoff_Callback}, 'BackgroundColor', 'white');
    rndmaxbutton = uicontrol('Style', 'edit', 'String','1.0', 'Position',...
        [149,163,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    rndlabel = uicontrol('Style', 'text', 'String', 'Roundness ', 'FontSize',9, ...
        'Position', [15,163,68,25]); 
    
    circminbutton = uicontrol('Style', 'edit', 'String','0.7', 'Position',...
        [91,127,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    circmaxbutton = uicontrol('Style', 'edit', 'String','1.0', 'Position',...
        [149,127,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    circlabel = uicontrol('Style', 'text', 'String', 'Circularity ', 'FontSize',9, ...
        'Position', [16,126,67,25]);
    
    diamminbutton = uicontrol('Style', 'edit', 'String','0', 'Position',...
        [91,235,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    diammaxbutton = uicontrol('Style', 'edit', 'String','100', 'Position',...
        [149,235,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    diamlabel = uicontrol('Style', 'text', 'String', 'Diameter ', 'FontSize',9, ...
        'Position', [31,236,52,23]);
    
    feretminbutton = uicontrol('Style', 'edit', 'String','0', 'Position',...
        [91,199,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    feretmaxbutton = uicontrol('Style', 'edit', 'String','100', 'Position',...
        [149,199,42,24], 'Callback', {@rdcutoff_Callback},'BackgroundColor', 'white');
    feretlabel = uicontrol('Style', 'text', 'String', 'Feret ', 'FontSize',9, ...
        'Position', [31,200,52,24]);
    
    %BUTTONS AT THE TOP 
    %select files button 
    selectfilebutton = uicontrol('Style', 'pushbutton', 'String','Select Folder','FontSize',8,...
        'Position',[10,351,78,23],'Callback', {@selectfile_Callback}); 
    %[left, bottom, width, height]
    
   %save figure
    svfigbutton = uicontrol('Style', 'pushbutton', 'String','Save .fig','FontSize',8,...
        'Position',[96,351,63,23],'Callback', {@svfig_Callback}); 
    
    %save jpg 
    svimagebutton = uicontrol('Style', 'pushbutton', 'String','Save .tiff','FontSize',8,...
        'Position',[166,351,63,23],'Callback', {@svimage_Callback}); 
    
    %sample ID
    sampleID = uicontrol('Style', 'edit', 'String','Sample Name', 'Position',...
        [91,322,99,22], 'Callback', {@sampleID_Callback},'BackgroundColor', 'white');
    sampleIDlabel = uicontrol('Style', 'text', 'String', 'Sample ID', 'FontSize',9, ...
        'Position', [16,317,68,25]); 
    
    %use compliment 
    checkboxcomp = uicontrol('Style', 'checkbox', 'String', 'black background','Position',[80,295,110,25],'Callback', {@donothing_Callback} );
  
    set(0, 'currentfigure', f3);
    hf3a3 = axes(); 
    imagelistlabel = uicontrol('Style', 'text', 'String', 'Image:', 'FontSize',9, ...
        'Position', [12,17,77,24]); 
    imagelistpop = uicontrol('Style','popupmenu', 'String','imagelist','Position',[80 19 93 22], 'Callback', {@imagelistpop_Callback});
    
    oldthresholdlabel = uicontrol('Style', 'text', 'String', 'Current Threshold:', 'FontSize',9, ...
        'Position', [183,17,107,19]); 
    oldthreshold = uicontrol('Style', 'text', 'String', '0.0', 'FontSize',9, ...
        'Position', [297, 23, 43, 12]); 
    
    newthresholdlabel = uicontrol('Style', 'text', 'String', 'New Threshold', 'FontSize',9, ...
        'Position', [350, 23, 99, 14]); 
    changethreshold = uicontrol('Style', 'edit', 'String', '0.0', 'FontSize',9, ...
        'Position', [455,17,53,21], 'Callback', {@changethreshold_Callback}); 
  
%%%%%%%%%GUI Initilization %%%%%%%%%%%%%%%
handles = get(f1, 'Children'); 
set(handles,'Units', 'normalized'); 
set(f1, 'Visible', 'on'); 
set(f1, 'Name', 'particle sizing')



%%%%%%%%Global Variable Definitions%%%%%%%%%%%%%%%
sourceDir = [];  
feretTotal = []; 
roundTotal = []; 
circTotal = []; 
diameter = []; 
feretTotalEdited = []; 
roundTotalEdited = []; 
diametersEdited = []; 
circTotalEdited = []; 
anno1 = [];
anno2 = []; 
anno3 = []; 
fit1plot = []; 
fit2plot = []; 
dtimages = []; 


%%%%%%%%%%Callbacks%%%%%%%%%%%%%%%%%%%

%general structure should have an 'update' function if any of the field are
%changes it updates the edited data by checking whether the original data
%fits with all the different filters using a isinbounds function
%the fitting boxes is only relevant to the plotting functions

%load the data after after clicking on get folder button 
    function  selectfile_Callback(source, eventdata)
        % Edit by Xingzhi Wang 07/16/2019 to enable saving previously
        % selected directory
        
        if isfile('selected_dir.mat')
            load('selected_dir.mat')
        else
            init_dir = path;
        end
        sourceDir = uigetdir(init_dir, 'select source folder'); 
        if sourceDir~=0
            init_dir = sourceDir;
            save('selected_dir.mat', 'init_dir')
        end       
        
        %sourceDir = uigetdir(path, 'select source folder'); 
        
        % End of Xingzhi Wang edit
        
        
        if get(checkboxcomp, 'Value'); 
            dtimages = fanalyzeimages(sourceDir, true);
        else 
            dtimages = fanalyzeimages(sourceDir, false);
        end 
        %the fanalyzeimages returns a stucture dtimage(j).field contains
        %various fields with information on the jth image 
        
        readImageStructure(); 
         
        
        popuplabel = []; 
        for h = 1:length(dtimages); 
            if (h < 10)
               popuplabel = [popuplabel; ['0', num2str(h)]];   
            else
            popuplabel = [popuplabel; num2str(h)]; 
            end
        end 
        set(imagelistpop, 'String', popuplabel); 
        
        updateGlobalData();
    end 
    
    function readImageStructure()
        diameter = []; 
        feretTotal = []; 
        roundTotal = []; 
        circTotal = []; 
        
        for j = 1:length(dtimages); 
                diameter = [diameter;[dtimages(j).properties(:).diameter]'];  
                feretTotal = [feretTotal;[dtimages(j).properties(:).nmFeret]']; 
                roundTotal = [roundTotal;[dtimages(j).properties(:).Round]']; 
                circTotal = [circTotal;[dtimages(j).properties(:).Circ]'];
        end
    end 

    function updateGlobalData()
        logarray = inboundaries(diameter, feretTotal, roundTotal, circTotal); 
        feretTotalEdited = feretTotal(logarray); 
        roundTotalEdited = roundTotal(logarray); 
        diametersEdited = diameter(logarray); 
        circTotalEdited = circTotal(logarray); 
        plot1(diametersEdited); 
        plot2(feretTotalEdited); 
        plot3(circTotal, roundTotal, diameter);
        plotimage(); 
    end 
    
    function [log] = inboundaries(diam, feret, round, circ)
        %pass a vector of data  diameter, feret, round, circularity 
        roundmax = str2num(get(rndmaxbutton, 'String')); 
        roundmin = str2num(get(rndminbutton, 'String')); 
        circmax = str2num(get(circmaxbutton, 'String')); 
        circmin = str2num(get(circminbutton, 'String')); 
        feretmax = str2num(get(feretmaxbutton, 'String')); 
        feretmin = str2num(get(feretminbutton, 'String')); 
        diammax = str2num(get(diammaxbutton, 'String')); 
        diammin = str2num(get(diamminbutton, 'String')); 
        
        log =  (diammin < diam) & ( diam < diammax) &...
               (feretmin < feret) & ( feret < feretmax) &...
               (circmin < circ) & ( circ < circmax) &...
               (roundmin < round) & ( round < roundmax); 
    end 
    
    function imagelistpop_Callback(source, eventdata)
        plotimage(); 
    end


    function plotimage()
        j = get(imagelistpop, 'Value'); 
        set(0, 'currentfigure', f3);
        axes(hf3a3)
        % whiteImage = 255 * dtimages(j).original;
        % imshow(whiteImage, 'Parent', hf3a3);
        imshow(dtimages(j).original, 'Parent', hf3a3); 
        hold on 
        angles = [];
        phi = [];
        greenAngle = [];
        cyanAngle = [];
        redAngle = [];
        outlier = [];
        % X1 = [];
        % X2 = [];
        % Y1 = [];
        % Y2 = [];

        t = linspace(0,2*pi,50);
        x1 = [];
        y1 = [];
        x2 = [];
        y2 = [];


        B = bwboundaries(dtimages(j).binary);   
                for k=1:(length(dtimages(j).properties))
                    if (inboundaries(dtimages(j).properties(k).diameter, dtimages(j).properties(k).nmFeret,...
                            dtimages(j).properties(k).Round, dtimages(j).properties(k).Circ)); 
                    b2 = B{k};

                    a = dtimages(j).properties(k).MajorAxisLength/2;
                    b = dtimages(j).properties(k).MinorAxisLength/2;
                    Xc = dtimages(j).properties(k).Centroid(1);
                    Yc = dtimages(j).properties(k).Centroid(2);
                    phi(k) = deg2rad(-dtimages(j).properties(k).Orientation);
                    
                    x = Xc + a*cos(t)*cos(phi(k)) - b*sin(t)*sin(phi(k));
                    y = Yc + a*cos(t)*sin(phi(k)) + b*sin(t)*cos(phi(k));

                    phi(k) = rad2deg(phi(k));
                    

                    if (phi(k) <= 30 && phi(k) >= -30) || (phi(k) >= 150 && phi(k) <= -150)
                        if (phi(k) >= 150 && phi(k) <= -150)
                            phi(k) = 180 + phi(k);
                        end
                        x1(k) = Xc - abs(a)*cosd(phi(k));
                        x2(k) = Xc + abs(a)*cosd(phi(k));
                        y1(k) = Yc - abs(a)*sind(phi(k));
                        y2(k) = Yc + abs(a)*sind(phi(k));
                        cyanAngle = [cyanAngle, phi(k)];

                        quiver(x1(k),y1(k),x2(k)-x1(k),y2(k)-y1(k),0, 'Color', 'c', 'LineWidth',3, 'MaxHeadSize', 10)


                    elseif (phi(k) >= 30 && phi(k) <= 90) || (phi(k) <= -90 && phi(k) >= -150)
                        if (phi(k) <= -90 && phi(k) >= -150)
                            phi(k) = 180 + phi(k);
                        end
                        x1(k) = Xc - abs(a)*cosd(phi(k));
                        x2(k) = Xc + abs(a)*cosd(phi(k));
                        y1(k) = Yc - abs(a)*sind(phi(k));
                        y2(k) = Yc + abs(a)*sind(phi(k));
                        greenAngle = [greenAngle, phi(k)];

                        quiver(x1(k),y1(k),x2(k)-x1(k),y2(k)-y1(k),0, 'Color', 'g', 'LineWidth',3, 'MaxHeadSize', 10)
                    elseif (phi(k) >= 90 && phi(k) <= 150) || (phi(k) <= -30 && phi(k) >= -90)
                        if (phi(k) <= -30 && phi(k) >= -90)
                            phi(k) = 180 + phi(k);
                        end
                        x1(k) = Xc - abs(a)*cosd(phi(k));
                        x2(k) = Xc + abs(a)*cosd(phi(k));
                        y1(k) = Yc - abs(a)*sind(phi(k));
                        y2(k) = Yc + abs(a)*sind(phi(k));

                        redAngle = [redAngle, phi(k)];
                        quiver(x1(k),y1(k),x2(k)-x1(k),y2(k)-y1(k),0, 'Color', 'r', 'LineWidth',3, 'MaxHeadSize', 10)
                    else
                        x1(k) = Xc - abs(a)*cosd(phi(k));
                        x2(k) = Xc + abs(a)*cosd(phi(k));
                        y1(k) = Yc - abs(a)*sind(phi(k));
                        y2(k) = Yc + abs(a)*sind(phi(k));
                        outlier = [outlier, phi(k)];
                        % quiver(x1(k),y1(k),x2(k)-x1(k),y2(k)-y1(k),0, 'Color', 'k', 'LineWidth',3, 'MaxHeadSize', 10)

                    end

                    %This below plots an ellipse and not the exact edge of the particle
					% plot(x,y,'m','Linewidth',0.5)

                    %This plot two points along the major axis for knowing the start and end points
					% plot(x1(k),y1(k), 'bo',x2(k),y2(k), 'go',x,y,'r','Linewidth',1)
                    
                    plot(hf3a3, b2(:,2),b2(:,1),'k','LineWidth',1);

                    else 
                    b2 = B{k};
					% plot(x,y,'b','Linewidth',0.5)
                    plot(hf3a3, b2(:,2),b2(:,1),'b','LineWidth',1);
                    end 
                    
                end
        
        %This is to find the area of the dots compared to the total area to get the packing density
        figure(4)
        imshow(dtimages(j).binary)
        % % stats = regionprops('table',dtimages(j).binary,'Area');
        % Dotstotalarea = bwarea(dtimages(j).binary)
        % Imagetotalarea = numel(dtimages(j).binary)

        % fraction = Dotstotalarea / Imagetotalarea


        title(hf3a3, dtimages(j).name); 
        set(oldthreshold, 'String', num2str(dtimages(j).threshold)); 
        
        %This is for the stats of the angle between the vectors
        MeanAngle = abs(mean(cyanAngle))
        MeanAnglegreen = abs(mean(greenAngle))
        MeanAnglered = abs(mean(redAngle))
        
        cyanAngle = cyanAngle - MeanAngle;
        greenAngle = greenAngle - MeanAngle;
        redAngle = redAngle - MeanAngle;
        

        figure(5)
        hold on
        h_cyan = histfit(cyanAngle);
        h_cyan(1).FaceColor = 'c';
        h_cyan(1).EdgeColor = 'k';
        h_cyan(1).EdgeAlpha = 0.5;
        h_cyan(2).Color = 'k';

        h_green = histfit(greenAngle);
        h_green(1).FaceColor = 'g';
        h_green(1).EdgeColor = 'k';
        h_green(1).EdgeAlpha = 0.5;
        h_green(2).Color = 'k';

        h_red = histfit(redAngle);
        h_red(1).FaceColor = 'r';
        h_red(1).EdgeColor = 'k';
        h_red(1).EdgeAlpha = 0.5;
        h_red(2).Color = 'k';

        ylabel('Counts', 'FontSize', 16)
        xlabel('\phi', 'FontSize', 16)
        set(gca, 'XLim', [-60 180])
        grid on
        grid minor
        ax = gca;
		ax.GridColor = 'k';
		ax.GridLineStyle = '--';
		ax.GridAlpha = 0.5;

    end 
        
    function plot1(data)
        cla(ha1)
        [g, x] = hist(data, 34); 
        bar(ha1, x, g/trapz(x,g), 1, 'FaceColor', 'w', 'EdgeColor', 'k');
        xlabel(ha1,'Diameter'); 
        ylabel(ha1,'Probability Density'); 
       
            %text(2, (max(pdfunction)-0.05), string);
            %probability distribution 
            checked = get(checkbox1, 'Value'); 
            if checked
                axes(ha1); 
                hold on 
                exclusion1 = diametersEdited<str2num(get(fit1max, 'String'));
                diametersExcluded = diametersEdited(exclusion1); 
                exclusion2 = diametersExcluded>str2num(get(fit1min, 'String'));
                diametersExcluded = diametersExcluded(exclusion2); 
                pd = fitdist(diametersExcluded, 'normal'); 
                mean1 = pd.mu; 
                std1 = pd.std; 
                [g, x] = hist(diametersEdited, 34);
                pdfunction = pdf(pd,x); 
                fit1plot = plot(ha1, x, pdfunction, 'r', 'LineWidth', 2); 
                hold off 
            end      
            if exist('anno1')
                delete(anno1);  
            end   
        label1 = ['d = ', num2str(mean(data), 4), ' +/- ', num2str(std(data), 3), ' nm '];
            if checked 
                labelfit1 = ['fit d = ', num2str(mean1, 4), '+/-', num2str(std1, 3), 'nm'];
                anno1 = annotation(f1, 'textbox',get(ha1,'Position'),'String', [label1,char(10),labelfit1],'LineStyle','none');
                legend(ha1, get(sampleID,'String'), 'fit');
            else
                anno1 = annotation(f1, 'textbox',get(ha1,'Position'),'String', label1,'LineStyle','none');
                legend(ha1, get(sampleID,'String'));
            end
          
        numberofparticles = length(data); 
        title(ha1, ['Estimated Diameters From Area with ', num2str(numberofparticles), ' particles']); 
    end 
    
    function plot2(data)
       cla(ha2)
       [g, x] = hist(data, 34); 
       bar(ha2, x, g/trapz(x,g), 1, 'FaceColor', 'w', 'EdgeColor', 'k'); 
       xlabel(ha2, 'Diameter'); 
       ylabel(ha2, 'Probability Density'); 
       label2 = ['d = ', num2str(mean(data), 4), ' +/- ', num2str(std(data), 3), ' nm'];
             %text(2, (max(pdfunction)-0.05), string);
             checked = get(checkbox2, 'Value'); 
             if checked
                axes(ha2); 
                hold on 
                exclusion1 = feretTotalEdited<str2num(get(fit2max, 'String'));
                feretExcluded = feretTotalEdited(exclusion1); 
                exclusion2 = feretExcluded>str2num(get(fit2min, 'String'));
                feretExcluded = feretExcluded(exclusion2); 
                pd = fitdist(feretExcluded, 'normal'); 
                mean2 = pd.mu; 
                std2 = pd.std; 
                [g, x] = hist(feretTotalEdited, 34);
                pdfunction = pdf(pd,x); 
                fit2plot = plot(ha2, x, pdfunction, 'r', 'LineWidth', 2); 
                hold off
             end 
       if exist('anno2')
            delete(anno2);  
        end  
       
        if checked 
                labelfit2 = ['fit d = ', num2str(mean2, 4), '+/-', num2str(std2, 3), 'nm'];
                anno2 = annotation(f1, 'textbox',get(ha2,'Position'),'String', [label2,char(10),labelfit2],'LineStyle','none');
                legend(ha2, get(sampleID,'String'),'fit'); 
            else
                anno2 = annotation(f1, 'textbox',get(ha2,'Position'),'String', label2,'LineStyle','none');
                legend(ha2, get(sampleID,'String'));
            end
        
       title(ha2, 'Feret Diameter'); 
    end

    function plot3(datax, datay, color) 
        cla(ha3); 
        %hist(ha3, data, 30); 
        %[g, x] = hist(data, 30); 
        %bar(ha3, x, g/trapz(x,g), 1); 
        %color manipulation
        
        scatter(ha3, datax, datay,[], color, 'SizeData', 4, 'MarkerFaceColor', 'flat'); 
        set(0, 'currentfigure', f1);
        axes(ha3); 
        colorbar; 
        xlabel(ha3, 'Circularity'); 
        ylabel(ha3, 'Roundness'); 
%             checked = get(checkbox3, 'Value'); 
%              if checked
%                 fprintf('checked evaluating \n')
%                 axes(ha3); 
%                 hold on 
%                 exclusion1 = data<str2num(get(fit3max, 'String'));
%                 roundExcluded = data(exclusion1); 
%                 exclusion2 = roundExcluded>str2num(get(fit3min, 'String'));
%                 roundExcluded = roundExcluded(exclusion2); 
%                 pd = fitdist(roundExcluded, 'normal'); 
%                 mean3 = pd.mu; 
%                 std3 = pd.std; 
%                 [g, x] = hist(data, 34);
%                 pdfunction = pdf(pd,x); 
%                 fit3plot = plot(ha3, x, pdfunction, 'r', 'LineWidth', 2); 
%                 hold off
%              end 
        if ~(isnan(anno3))
            delete(anno3)
        end
        %label3 = ['r = ', num2str(mean(data), 4), ' +/- ', num2str(std(data), 3)];
        %anno3 = annotation(f1, 'textbox',get(ha3,'Position'),'String', label3,'LineStyle','none');
        title(ha3, 'Morphological Data') 
    end


    function rdcutoff_Callback(source,eventdata)
        updateGlobalData(); 
    end 

    function checkbox1_Callback(source, eventdata)
    plot1(diametersEdited) 
    end
    function fit1min_Callback(source, eventdata)
    checkbox1_Callback(checkbox1)
    end
    function fit1max_Callback(source, eventdata)
    checkbox1_Callback(checkbox1)
    end

    function checkbox2_Callback(source, eventdata)
    plot2(feretTotalEdited)     
    end
    function fit2min_Callback(source, eventdata)
    checkbox2_Callback(checkbox2)
    end
    function fit2max_Callback(source, eventdata)
    checkbox2_Callback(checkbox2)
    end

    function sampleID_Callback(source, eventdata)
    plot1(diametersEdited); 
    plot2(feretTotalEdited);     
    end

    function svfig_Callback(source, eventdata)
    savefig(f1, fullfile(sourceDir, get(sampleID, 'String'))); 
    end

    function svimage_Callback(source, eventdata)
    set(f1,'PaperPositionMode','auto')
    print(f1,'-r600','-opengl', '-dtiff', fullfile(sourceDir, get(sampleID, 'String')))
    end

    function changethreshold_Callback(source, eventdata)
        j = get(imagelistpop, 'Value');
        manualthresh = str2num(get(changethreshold,'String')); 
        [BW, filtered, thresh] = fthreshold1(dtimages(j).edited,manualthresh,false);
        dtimages(j).binary = BW; 
        dtimages(j).threshold = thresh; 
        jthImageProperties = fCalculateProperties(BW, dtimages(j).magnification); 
        dtimages(j).properties = jthImageProperties;
        
        readImageStructure(); 
        updateGlobalData(); 
    end 

    function donothing_Callback(source, eventdata)
    end

end

