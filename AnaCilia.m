function AnaCilia()
% [Description]


%% == Parameter setting ==
clc;clearvars;
p=parameter_setting();
warning('off','all');
% -------------------------------------------------
%% == Import image file ==
folder = fileparts(which(p.FileName));
fullFileName = fullfile(folder,p.FileName);
if ~exist(fullFileName, 'file')
	if ~exist(baseFileName, 'file')
		return;
    end
	fullFileName = p.FileName;
end

%% == All tiff in image ==
[IMseries,~]=ImageTiffSeries(fullFileName);
max_IMi=size(IMseries,2)/4;[p.ny,p.nx]=size(IMseries{1,1});

if mod(size(IMseries,2),4)~=0
else
    for IMi=1:max_IMi

    %% == Analysis cell image ==
    ori_Image = imread(fullFileName,IMi); % ori_Image: NxM unit16

    binaryImage=PreImage(ori_Image,p,p.thresholdValue);
    [binaryImage,blob,blob2]=BlobImage(binaryImage,ori_Image,p);

    %% == Analysis cilia image ==
    cili_Image = imread(fullFileName,IMi+max_IMi*2); % ori_Image: NxM unit16
    eval(['figure(''visible'',''' p.displayimage ''');imagesc(cili_Image);']); 

    [binaryImage_c,cili_blob]=cilia(cili_Image,p);

    %% == Integrated cell and cilia analysis ==
    [blob2,cili_blob]=cell_cilia(blob2,cili_blob,p);
    [BW1,BW0,result]=draw_cell_cilia(binaryImage,blob2,cili_blob,cili_Image,ori_Image,binaryImage_c,p);

    eval(['save(''result_' num2str(IMi) '.mat'',''IMseries'',''blob'',''blob2'',''cili_blob'',''BW1'',''BW0'',''result'',''ori_Image'',''p'');']);
    end
end
end
function p=parameter_setting()
p.FileName='WT10%FBS60X1.tif';
p.rgb=3;
p.rgb2=23;
p.thresholdValue = 150;
p.cilia_thresholdValue = 1500;

p.se_tophat=30;
p.se_open=15;
p.se_erode=5;
p.num_erode=1;
p.maxgray_o=65536;
p.maxgray_new=4000;
p.threhold_step=100;

p.displayimage='off';   % 'off' or 'on'
p.displayimage1='on';


% programming process control
end
function binaryImage=PreImage(IM,p,thresh)

%% == 1. Image Top-hat filtering ==

originalImage=filter_tophat(IM,p.se_tophat,1);
%figure;imshow(originalImage);


%% == 2. Draw histogram and use threshold value to identify object ==
[~,binaryImage]=imhist_re(originalImage,p.maxgray_o,p.maxgray_new,p,thresh);


%% == 3. Morphologically open image ==
binaryImage=filter_open(binaryImage,p.se_open,1);
%figure;imshow(binaryImage);


%% == 4. Morphologically erode image ==
%binaryImage=filter_erode(binaryImage,p.se_erode,p.num_erode);
%figure;imshow(binaryImage);

    function [pixelC,binaryImage]=imhist_re(originalImage,maxgray_o,maxgray_new,p,thresh)
        %Draw histogram figure and use threshold value to cut image
        
        [pix, gra] = imhist(originalImage,maxgray_o);
        pixelC=pix(1:maxgray_new,1);grayl=gra(1:maxgray_new,1);
        %figure;bar(pixelC);
        %hold on;
        %maxYValue = ylim;
        %line([p.thresholdValue, p.thresholdValue], maxYValue, 'Color', 'r');
        binaryImage = originalImage > thresh;
        %hold off;
    end
end
function [binaryImage2,blob,blob2]=BlobImage(binaryImage,originalImage,p)
% Image process for one image

%% == 5. Get all the blob properties and draw boundary ==
blob=BlobInfo(binaryImage,originalImage,p);
%binaryImage2=Seperate_blob(blob,binaryImage,p);
binaryImage2=Seperate_blob(blob,p);
blob2=BlobInfo(binaryImage2,originalImage,p);
%BlobBounrady(binaryImage2,originalImage,p);

    function [blob,blobMeasurements]=BlobInfo(BW,IM,p)
        % Get all the blob properties and draw boundary
        % BW: binaryImage
        % IM: originalImage
        
        labeledImage = bwlabel(BW, 4);
        blobMeasurements = regionprops(labeledImage, IM, 'all');
        numberOfBlobs = size(blobMeasurements, 1);
        
        eval(['figure(''visible'',''' p.displayimage ''');imshow(BW);']); 
        %figure('visible','off');imshow(BW); 
        hold on;
        for k = 1 : numberOfBlobs
            blob{k}.thisBlobsPixels = blobMeasurements(k).PixelIdxList;
            blob{k}.meanGL = mean(IM(blob{k}.thisBlobsPixels));	
            blob{k}.blobArea = blobMeasurements(k).Area;
            blob{k}.blobPerimeter = blobMeasurements(k).Perimeter;
            blob{k}.blobCentroid = blobMeasurements(k).Centroid;
            blob{k}.Image = blobMeasurements(k).Image;
            blob{k}.BoundingBox = blobMeasurements(k).BoundingBox;
            blob{k}.PixelList = blobMeasurements(k).PixelList;
            blob{k}.PixelValues = blobMeasurements(k).PixelValues;
            textFontSize = 8;labelShiftX = -2;
            text(blob{k}.blobCentroid(1) + labelShiftX, blob{k}.blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','r');
        end
        hold off;
    end
    function BlobBounrady(BW,IM,p)
            %figure;imshow(IM);
            eval(['figure(''visible'',''' p.displayimage ''');imshow(IM);']); 
            hold on;
            boundaries = bwboundaries(BW);
            numberOfBoundaries = size(boundaries, 1);
            for k = 1 : numberOfBoundaries
                thisBoundary = boundaries{k};
                plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);  
            end
            hold off;
    end
end
function [binaryImage_c,blob]=cilia(IM,p)

%% == 1. Draw histogram and use threshold value to identify object ==
[~,binaryImage_c]=imhist_re(IM,p.maxgray_o,p.maxgray_new);
%figure;imshow(binaryImage);

%% == 2. Get all the blob properties and draw boundary ==
blob=BlobInfo(binaryImage_c,IM,p);

    function [pixelC,binaryImage]=imhist_re(originalImage,maxgray_o,maxgray_new)
        %Draw histogram figure and use threshold value to cut image
        
        [pix, gra] = imhist(originalImage,maxgray_o);
        pixelC=pix(1:maxgray_new,1);grayl=gra(1:maxgray_new,1);
        %figure;bar(pixelC);
        eval(['figure(''visible'',''' p.displayimage ''');bar(pixelC);']); 
        hold on;
        maxYValue = ylim;
        line([p.cilia_thresholdValue, p.cilia_thresholdValue], maxYValue, 'Color', 'r');
        binaryImage = originalImage > p.cilia_thresholdValue;
        hold off;
    end
    function [blob,blobMeasurements]=BlobInfo(BW,IM,p)
        % Get all the blob properties and draw boundary
        % BW: binaryImage
        % IM: originalImage
        
        labeledImage = bwlabel(BW, 4);
        blobMeasurements = regionprops(labeledImage, IM, 'all');
        numberOfBlobs = size(blobMeasurements, 1);
        
        %figure;imshow(BW);
        eval(['figure(''visible'',''' p.displayimage ''');imshow(BW);']); 
        hold on;
        for k = 1 : numberOfBlobs
            blob{k}.thisBlobsPixels = blobMeasurements(k).PixelIdxList;
            blob{k}.meanGL = mean(IM(blob{k}.thisBlobsPixels));	
            blob{k}.blobArea = blobMeasurements(k).Area;
            blob{k}.blobPerimeter = blobMeasurements(k).Perimeter;
            blob{k}.blobCentroid = blobMeasurements(k).Centroid;
            blob{k}.Image = blobMeasurements(k).Image;
            blob{k}.BoundingBox = blobMeasurements(k).BoundingBox;
            blob{k}.PixelList = blobMeasurements(k).PixelList;
            blob{k}.PixelValues = blobMeasurements(k).PixelValues;
            textFontSize = 8;labelShiftX = -2;
            text(blob{k}.blobCentroid(1) + labelShiftX, blob{k}.blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','r');
        end
        hold off;
    end
end
function [blob,cili_blob]=cell_cilia(blob,cili_blob,p)

nc_max=size(blob,2); % number of cell number
ni_max=size(cili_blob,2); % number of cilia number
c_center=zeros(nc_max,3);i_center=zeros(ni_max,3);

% create center matrix for cell and cilia
for nc=1:nc_max
    c_center(nc,1:2)=blob{1,nc}.blobCentroid;
    blob{1,nc}.ciliaID=[];
end
for ni=1:ni_max
    i_center(ni,1:2)=cili_blob{1,ni}.blobCentroid;
end
% compare center matrixs
D = pdist2(c_center(:,1:2),i_center(:,1:2));
[d1,nd1]=min(D);i_center(:,3)=nd1';
%
for ni=1:ni_max
    c_center(i_center(ni,3),3)=c_center(i_center(ni,3),3)+1;
    cili_blob{1,ni}.cellID=i_center(ni,3);
    blob{1,i_center(ni,3)}.ciliaID=[blob{1,i_center(ni,3)}.ciliaID ni];  
end
for nc=1:nc_max
    blob{1,nc}.ciliaN=c_center(nc,3);
end

end
function [BW1,BW0,result]=draw_cell_cilia(binaryImage,blob,cili_blob,cili_Image,ori_Image,binaryImage_c,p)
    % BW1: cell with cilia, BW0: cell without cilia
    BW1=im2bw(zeros(size(binaryImage)));BW0=binaryImage;n_cell1=0;
    for k = 1 : size(blob,2)
        if blob{1,k}.ciliaN>0
            pixelList_image=[blob{1,k}.PixelList(:,1) blob{1,k}.PixelList(:,2)];
            for kk=1:size(pixelList_image,1)
                BW1(pixelList_image(kk,2),pixelList_image(kk,1))=1;
                BW0(pixelList_image(kk,2),pixelList_image(kk,1))=0;
            end
            n_cell1=n_cell1+1;
        end
    end
    table1=zeros(n_cell1,3);
    %figure;imshow(BW1);
    eval(['figure(''visible'',''' p.displayimage ''');imshow(BW1);']); 
    textFontSize = 8;labelShiftX = -2;
    hold on;k0=0;
    for k = 1 : size(blob,2)
        if blob{1,k}.ciliaN>0
            k0=k0+1;
            text(blob{k}.blobCentroid(1) + labelShiftX, blob{k}.blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','r');
            table1(k0,1)=k;
            temp=0;n_temp=0;
            for kk=1:size(blob{k}.ciliaID,2)
                pi_position=cili_blob{1,blob{k}.ciliaID(1,kk)}.PixelList;
                for kkk=1:size(pi_position,1)
                    n_temp=n_temp+1;temp=temp+cili_Image(pi_position(kkk,2),pi_position(kkk,1));
                end
            end
            table1(k0,2)=temp/n_temp;table1(k0,3)=n_temp;
        end
    end
    hold off;
    %figure;imshow(BW0);
    eval(['figure(''visible'',''' p.displayimage ''');imshow(BW0);']); 
    % cell-cilia result
    result.n_cell1=n_cell1;
    result.n_cell0=size(blob,2)-n_cell1;
    result.pect_cili=n_cell1/size(blob,2);
    result.table1=table1;
    
    % Draw boundary-marker-cell with cilia figure
    %figure;imagesc(ori_Image);
    eval(['figure(''visible'',''' p.displayimage1 ''');imagesc(ori_Image);']); 
    colormap parula;
    hold on;
    for k = 1 : size(table1,1)
        textFontSize = 8;labelShiftX = -2;
        text(blob{1,table1(k,1)}.blobCentroid(1) + labelShiftX, blob{1,table1(k,1)}.blobCentroid(2), num2str(table1(k,1)), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','r');
    end
    hold on;
    boundaries = bwboundaries(BW1);
    numberOfBoundaries = size(boundaries, 1);
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);  
    end
    hold on;
    boundaries0 = bwboundaries(binaryImage_c);
    numberOfBoundaries0 = size(boundaries0, 1);
    for k = 1 : numberOfBoundaries0
        thisBoundary0 = boundaries0{k};
        plot(thisBoundary0(:,2), thisBoundary0(:,1), 'w', 'LineWidth', 3);  % Cilia color
    end
    title('Cell with cilia');
    hold off;
    
    % Draw all boundary-marker-cell figure
    %figure;imagesc(ori_Image);
    eval(['figure(''visible'',''' p.displayimage1 ''');imagesc(ori_Image);']); 
    colormap parula;
    hold on;
    for k = 1 : size(blob,2)
        textFontSize = 8;labelShiftX = -2;
        text(blob{1,k}.blobCentroid(1) + labelShiftX, blob{1,k}.blobCentroid(2), num2str(k), 'FontSize', textFontSize, 'FontWeight', 'Bold','Color','r');
    end
    hold on;
    boundaries = bwboundaries(binaryImage);
    numberOfBoundaries = size(boundaries, 1);
    for k = 1 : numberOfBoundaries
        thisBoundary = boundaries{k};
        plot(thisBoundary(:,2), thisBoundary(:,1), 'g', 'LineWidth', 2);  
    end
    title('Total Cell');
    hold off;
    
    % Draw cilia bar figure
    %figure;bar(result.table1(:,2));
    eval(['figure(''visible'',''' p.displayimage ''');bar(result.table1(:,2));']); 
    x_label=[];
    for xl=1:size(result.table1,1)
        b=result.table1(xl,1);
        if size(num2str(100),2)==1
            bb=['  ' b];
        elseif size(num2str(100),2)==2
            bb=[' ' b];
        elseif size(num2str(100),2)==3
            bb=b;
        else
            bb='000';
        end
        x_label=[x_label;bb];
    end
    set(gca,'XTick',1:size(result.table1,1));set(gca,'XTickLabel',x_label);
    xlabel('Cell ID');ylabel('The cilia staining intensity');
    eval(['title(''#Cell with cilia = ' num2str(result.n_cell1) '   #Total Cell = ' num2str(result.n_cell1+result.n_cell0) '   Cilia cell ratio% = ' num2str(result.pect_cili*100) ''')']);
end

% filtering
function IM2=filter_erode(IM1,se_erode,n_erode)
% [ morphological erode filtering ]
% IM1,IM2: gray or binary image
% se_erode: disk size of erode areucturing element
% n_erode: repeat number of erode filtering
    if n_erode==0
    else
    se=strel('disk',se_erode);
        for n=1:n_erode
            IM1=imerode(IM1,se);
        end
    end
    IM2=IM1;
end
function IM2=filter_dilate(IM1,se_dilate,n_dilate)
% [morphological dilate filtering ]
% IM1,IM2: gray or binary image
% se_erode: disk size of dilate areucturing element
% n_erode: repeat number of dilate filtering

    if n_dilate==0
    else
    se=strel('disk',se_dilate);
        for n=1:n_dilate
            IM1=imdilate(IM1,se);
        end
    end
    IM2=IM1;
end
function IM2=filter_tophat(IM1,se_tophat,n_tophat)
% [ morphological top-hat filtering ]
% IM1,IM2: gray or binary image
% se_tophat: disk size of tophat areucturing element
% n_tophat: repeat number of tophat filtering

    if n_tophat==0
    else
    se=strel('disk',se_tophat);
        for i=1:n_tophat
            IM1=imtophat(IM1,se);
        end
    end
    IM2=IM1;
end
function IM2=filter_bothat(IM1,se_bothat,n_bothat)
% [ morphological bottom-hat filtering ]
% IM1,IM2: gray or binary image
% se_tophat: disk size of bothat areucturing element
% n_tophat: repeat number of bothat filtering

    if n_bothat==0
    else
    se=strel('disk',se_bothat);
        for i=1:n_bothat
            IM1=imbothat(IM1,se);
        end
    end
    IM2=IM1;
end
function IM2=filter_open(IM1,se_open,n_open)
% [ morphological open filtering ]
% IM1,IM2: gray or binary image
% se_open: disk size of open areucturing element
% n_open: repeat number of open filtering

    if n_open==0
    else
    se=strel('disk',se_open);
        for i=1:n_open
            IM1=imopen(IM1,se);
        end
    end
    IM2=IM1;
end
function IM2=filter_close(IM1,se_close,n_close)
% [ morphological close filtering ]
% IM1,IM2: gray or binary image
% se_open: disk size of close areucturing element
% n_open: repeat number of close filtering

    if n_close==0
    else
    se=strel('disk',se_close);
        for i=1:n_close
            IM1=imclose(IM1,se);
        end
    end
    IM2=IM1;
end

% Seperate blob
function binaryImage2=Seperate_blob(blob,p)
% [ image watershed transform ]
% binaryImage2=Seperate_blob(blob,binaryImage,p)

thed=0.9;
thed2=10;
mac_whilen=20;

binaryImage2=zeros(p.ny,p.nx);
for bi=1:size(blob,2);
    obj=0;threshold_o=p.thresholdValue;whilen=0;
    
    % Adaptive blob modification
    while obj==0&&whilen<mac_whilen
        
        B=blob{bi};
        if exist('BW_ada','var')==1
            IM2d=BW_ada;
        else
            IM2d=B.Image;
            IM3d=zeros(size(IM2d));
            for k=1:B.blobArea
                IM3d(B.PixelList(k,2)-B.BoundingBox(1,2)+0.5,B.PixelList(k,1)-B.BoundingBox(1,1)+0.5)=B.PixelValues(k,1);
            end
            pre_area=B.blobArea;
        end
        %figure;imagesc(IM2d);
        [BW2,BWedge]=edge_add(IM2d);
        k2=conv_point(BWedge,thed,IM2d);
        IM2d_blob=link_neck(IM2d,k2,BWedge,thed2);
        %figure;imagesc(IM2d_blob);
        labeledImage = bwlabel(IM2d_blob, 4);
        blobMeasurements = regionprops(labeledImage, IM3d, 'all');
        sum_solod=0;num_solod=0;tol_area=0;
        %IM3d=zeros(size(IM2d));
        for bb=1:size(blobMeasurements,1)
            tol_area=tol_area+blobMeasurements(bb).Area;
            sum_solod=sum_solod+blobMeasurements(bb).Solidity;
            num_solod=num_solod+1;           
            for k=1:blobMeasurements(bb).Area
                %IM3d(blobMeasurements(bb).PixelList(k,2),blobMeasurements(bb).PixelList(k,1))=blobMeasurements(bb).PixelValues(k,1);
            end
        end
        %figure;imagesc(IM3d);
        ave_solod=sum_solod/num_solod;
        if abs(tol_area-pre_area)/pre_area>0.3||tol_area==0
            obj=1;IM2d_blob=pre_IM2d_blob;IM3d=pre_IM3d;
        elseif ave_solod>0.95
            obj=1;
        else
            threshold_o=threshold_o+p.threhold_step;
            BW_ada=PreImage(IM3d,p,threshold_o);
        end
        whilen=whilen+1;pre_area=tol_area;
        pre_IM2d_blob=IM2d_blob;pre_IM3d=IM3d;
    end
   
    % Draw inaryImage2
    for kx=1:size(IM2d_blob,2)
        for ky=1:size(IM2d_blob,1)
            if IM2d_blob(ky,kx)==1
                binaryImage2(ky+B.BoundingBox(1,2)-0.5,kx+B.BoundingBox(1,1)-0.5)=1;
            end
        end
    end
    clear BW_ada;
end    
   binaryImage2=im2bw(binaryImage2); 
end
function output=water_fix(IM3d,L,i,j,r,threh,n,m)
    region_t2=IM3d(max(i-r,1):min(i+r,n),max(j-r,1):min(j+r,m));
    region_t1=L(max(i-r,1):min(i+r,n),max(j-r,1):min(j+r,m));
    region_t=[region_t1(:)';region_t2(:)'];region_t=(sortrows(region_t'))';
    tabL=tabulate(region_t(1,:));
    tabl_out=[];ii=0;
    for k=1:size(tabL,1)
        if tabL(k,2)~=0
            tabl_out=[tabl_out;[tabL(k,1) mean(mean(region_t(2,ii+1:ii+tabL(k,2))))]];
            ii=ii+tabL(k,2);
        end
    end
    if tabl_out(1,1)==0
        tabl_out=tabl_out(2:end,:);
    end
    if std(tabl_out(:,2))<threh
        [a,b]=min(abs(tabl_out(:,2)-IM3d(i,j)));
        output=tabl_out(b,1);
    else
        output=0;
    end
end
function [BW2,BWedge]=edge_add(BW1)
% Obtain edge x y coordinates from BW image
BW2=BW1;
BWedge0 = bwboundaries(BW1);BWedge=[BWedge0{1,1}(:,2) BWedge0{1,1}(:,1)];
%{ 
Add edge in image boundaries
BW2=edge(BW1);
for i=1:size(BW1,1)
    for j=1:size(BW1,2)     
        if i==1||i==size(BW1,1)||j==1||j==size(BW1,2)
            if BW1(i,j)==1
                BW2(i,j)=1;
            end
        end     
    end
end
[xy(:,2),xy(:,1)]=find(BW2);n=size(xy,1); %find the point value=1 in BW2
xy=[xy zeros(n,1)];
BWedge=zeros(n,2);BWedge(1,:)=xy(1,1:2);xy(1,3)=1;
for k=2:n
    xy=xy(xy(:,3)~=1,:);
    A=distance_fix(xy(:,1:2),BWedge(k-1,:));
    [a,b]=min(A);xy(b,3)=1;
    BWedge(k,:)=xy(b,1:2);
end

BWedge=BWedge(BWedge(:,1)~=0,:);
%}
    function A=distance_fix(M,B_point)
        A=zeros(size(M,1),1);
        for k1=1:size(M,1)
            A(k1,1)=sqrt((M(k1,1)-B_point(1))^2+(M(k1,2)-B_point(2))^2);
        end  
    end
end
function k2=conv_point(BWedge,thed,IM2d)
% Obtain neck points k2 from BWedge xy coordinates with thed distance threshold

    % convex hull of BWedge: k
    k=convhull(BWedge);k2=[];
    if sum(diff(k)./abs(diff(k)))<0 % check order of k equal to BWedge
        k=k(end:-1:1);
    end
    
    %figure;plot(BWedge(:,1),BWedge(:,2));hold on;
    %plot(BWedge(k,1),BWedge(k,2));hold off;set(gca,'YDir','Reverse');
    
    % find neck points on cell boundary: k2
    k2=[];radius_bob=min([size(IM2d,1) size(IM2d,2)])/20;
    for i=2:size(k,1) % for each convex region
        
        xh=linspace(BWedge(k(i-1),1),BWedge(k(i),1),500); % linspace convex hull xy into 50 pieces: (xh,yh) 
        yh=linspace(BWedge(k(i-1),2),BWedge(k(i),2),500);
        xw=BWedge(k(i-1):k(i),1);yw=BWedge(k(i-1):k(i),2);w_list=k(i-1):k(i); % wall xy: (xw,yw)
        conv_d=zeros(1,size(xw,1));
        for j=1:size(xw,1)
            conv_d(1,j)=min(distance_fix([xh;yh]',[xw(j,1) yw(j,1)])); % find distance to hull for each wall point            
        end
        % Peak finder
        varargout = peakfinder(conv_d,(max(conv_d)-min(conv_d))/100,max(conv_d)/2,1); % find the local min points on Convex_dist
        [~,newk_n]=find(varargout~=1&varargout~=size(conv_d,2));
        if isempty(newk_n)~=1&&max(conv_d)>radius_bob
            newk_xy=varargout(newk_n);nn=size(newk_xy,2);
            for k3=1:nn
                if k(i-1)<k(i)  % specific region DOCH
                    newkk=newk_xy(k3)+k(i-1)-1;
                    %[~,newkk]=min((BWedge(newk_xy(k3)+k(i-1)-1,1)-x).^2+(BWedge(newk_xy(k3)+k(i-1)-1,2)-y).^2);
                    newkkk=[newkk;conv_d(newk_xy(k3))];
                elseif (k(i-1)>k(i))&&(k(i-1)-k(i)<=0.5*size(BWedge,1))  % specific region DOCH(inverse direction)
                    newkk=newk_xy(k3)-k(i-1)+1;
                    %[~,newkk]=min((cell.x_cell(start_convd-newk_xy(i)+1)-x).^2+(cell.y_cell(start_convd-newk_xy(i)+1)-y).^2);
                    newkkk=[newkk;conv_d(newk_xy(k3))];
                else  % specific region DOCH across(0,0)
                    % ===== below not check yet!!! === 
                    if newk_xy(i)<=(size(cell.x_cell,2)-start_convd)
                        [~,newkk]=min((cell.x_cell(newk_xy(i)+start_convd-1)-x).^2+(cell.y_cell(newk_xy(i)+start_convd-1)-y).^2);
                        newkkk=[newkk;conv_d(newk_xy(k3))];
                    else
                        [~,newkk]=min((cell.x_cell(newk_xy(i)-size(cell.x_cell,2)+start_convd)-x).^2+(cell.y_cell(newk_xy(i)-size(cell.x_cell,2)+start_convd)-y).^2);
                        newkkk=[newkk;conv_d(newk_xy(k3))];
                    end
                    % ================================
                end
                k2=[k2 newkkk];
            end
        end
        %[a,b]=max(conv_d); % max distance a to hull in each convex region, neck point: b
        %if a>thed
        %    k2=[k2 w_list(b)]; % k2: neck points
        %end
    end
    k2=shortk2(k2);
    function varargout = peakfinder(x0, sel, thresh, extrema)
        %PEAKFINDER Noise tolerant fast peak finding algorithm
        %   INPUTS:
        %       x0 - A real vector from the maxima will be found (required)
        %       sel - The amount above surrounding data for a peak to be
        %           identified (default = (max(x0)-min(x0))/4). Larger values mean
        %           the algorithm is more selective in finding peaks.
        %       thresh - A threshold value which peaks must be larger than to be
        %           maxima or smaller than to be minima.
        %       extrema - 1 if maxima are desired, -1 if minima are desired
        %           (default = maxima, 1)
        %   OUTPUTS:
        %       peakLoc - The indicies of the identified peaks in x0
        %       peakMag - The magnitude of the identified peaks
        %
        %   [peakLoc] = peakfinder(x0) returns the indicies of local maxima that
        %       are at least 1/4 the range of the data above surrounding data.
        %
        %   [peakLoc] = peakfinder(x0,sel) returns the indicies of local maxima
        %       that are at least sel above surrounding data.
        %
        %   [peakLoc] = peakfinder(x0,sel,thresh) returns the indicies of local 
        %       maxima that are at least sel above surrounding data and larger
        %       (smaller) than thresh if you are finding maxima (minima).
        %
        %   [peakLoc] = peakfinder(x0,sel,thresh,extrema) returns the maxima of the
        %       data if extrema > 0 and the minima of the data if extrema < 0
        %
        %   [peakLoc, peakMag] = peakfinder(x0,...) returns the indicies of the
        %       local maxima as well as the magnitudes of those maxima
        %
        %   If called with no output the identified maxima will be plotted along
        %       with the input data.
        %
        %   Note: If repeated values are found the first is identified as the peak
        %
        % Ex:
        % t = 0:.0001:10;
        % x = 12*sin(10*2*pi*t)-3*sin(.1*2*pi*t)+randn(1,numel(t));
        % x(1250:1255) = max(x);
        % peakfinder(x)
        %
        % Copyright Nathanael C. Yoder 2011 (nyoder@gmail.com)

        % Perform error checking and set defaults if not passed in
        error(nargchk(1,4,nargin,'struct'));
        error(nargoutchk(0,2,nargout,'struct'));

        s = size(x0);
        flipData =  s(1) < s(2);
        len0 = numel(x0);
        if len0 ~= s(1) && len0 ~= s(2)
            error('PEAKFINDER:Input','The input data must be a vector')
        elseif isempty(x0)
            varargout = {[],[]};
            return;
        end
        if ~isreal(x0)
            warning('PEAKFINDER:NotReal','Absolute value of data will be used')
            x0 = abs(x0);
        end

        if nargin < 2 || isempty(sel)
            sel = (max(x0)-min(x0))/4;
        elseif ~isnumeric(sel) || ~isreal(sel)
            sel = (max(x0)-min(x0))/4;
            warning('PEAKFINDER:InvalidSel',...
                'The selectivity must be a real scalar.  A selectivity of %.4g will be used',sel)
        elseif numel(sel) > 1
            warning('PEAKFINDER:InvalidSel',...
                'The selectivity must be a scalar.  The first selectivity value in the vector will be used.')
            sel = sel(1);
        end

        if nargin < 3 || isempty(thresh)
            thresh = [];
        elseif ~isnumeric(thresh) || ~isreal(thresh)
            thresh = [];
            warning('PEAKFINDER:InvalidThreshold',...
                'The threshold must be a real scalar. No threshold will be used.')
        elseif numel(thresh) > 1
            thresh = thresh(1);
            warning('PEAKFINDER:InvalidThreshold',...
                'The threshold must be a scalar.  The first threshold value in the vector will be used.')
        end

        if nargin < 4 || isempty(extrema)
            extrema = 1;
        else
            extrema = sign(extrema(1)); % Should only be 1 or -1 but make sure
            if extrema == 0
                error('PEAKFINDER:ZeroMaxima','Either 1 (for maxima) or -1 (for minima) must be input for extrema');
            end
        end

        x0 = extrema*x0(:); % Make it so we are finding maxima regardless
        thresh = thresh*extrema; % Adjust threshold according to extrema.
        dx0 = diff(x0); % Find derivative
        dx0(dx0 == 0) = -eps; % This is so we find the first of repeated values
        ind = find(dx0(1:end-1).*dx0(2:end) < 0)+1; % Find where the derivative changes sign

        % Include endpoints in potential peaks and valleys
        x = [x0(1);x0(ind);x0(end)];
        ind = [1;ind;len0];

        % x only has the peaks, valleys, and endpoints
        len = numel(x);
        minMag = min(x);


        if len > 2 % Function with peaks and valleys

            % Set initial parameters for loop
            tempMag = minMag;
            foundPeak = false;
            leftMin = minMag;

            % Deal with first point a little differently since tacked it on
            % Calculate the sign of the derivative since we taked the first point
            %  on it does not neccessarily alternate like the rest.
            signDx = sign(diff(x(1:3)));
            if signDx(1) <= 0 % The first point is larger or equal to the second
                ii = 0;
                if signDx(1) == signDx(2) % Want alternating signs
                    x(2) = [];
                    ind(2) = [];
                    len = len-1;
                end
            else % First point is smaller than the second
                ii = 1;
                if signDx(1) == signDx(2) % Want alternating signs
                    x(1) = [];
                    ind(1) = [];
                    len = len-1;
                end
            end

            % Preallocate max number of maxima
            maxPeaks = ceil(len/2);
            peakLoc = zeros(maxPeaks,1);
            peakMag = zeros(maxPeaks,1);
            cInd = 1;
            % Loop through extrema which should be peaks and then valleys
            while ii < len
                ii = ii+1; % This is a peak
                % Reset peak finding if we had a peak and the next peak is bigger
                %   than the last or the left min was small enough to reset.
                if foundPeak
                    tempMag = minMag;
                    foundPeak = false;
                end

                % Make sure we don't iterate past the length of our vector
                if ii == len
                    break; % We assign the last point differently out of the loop
                end

                % Found new peak that was lager than temp mag and selectivity larger
                %   than the minimum to its left.
                if x(ii) > tempMag && x(ii) > leftMin + sel
                    tempLoc = ii;
                    tempMag = x(ii);
                end

                ii = ii+1; % Move onto the valley
                % Come down at least sel from peak
                if ~foundPeak && tempMag > sel + x(ii)
                    foundPeak = true; % We have found a peak
                    leftMin = x(ii);
                    peakLoc(cInd) = tempLoc; % Add peak to index
                    peakMag(cInd) = tempMag;
                    cInd = cInd+1;
                elseif x(ii) < leftMin % New left minima
                    leftMin = x(ii);
                end
            end

            % Check end point
            if x(end) > tempMag && x(end) > leftMin + sel
                peakLoc(cInd) = len;
                peakMag(cInd) = x(end);
                cInd = cInd + 1;
            elseif ~foundPeak && tempMag > minMag % Check if we still need to add the last point
                peakLoc(cInd) = tempLoc;
                peakMag(cInd) = tempMag;
                cInd = cInd + 1;
            end

            % Create output
            peakInds = ind(peakLoc(1:cInd-1));
            peakMags = peakMag(1:cInd-1);
        else % This is a monotone function where an endpoint is the only peak
            [peakMags,xInd] = max(x);
            if peakMags > minMag + sel
                peakInds = ind(xInd);
            else
                peakMags = [];
                peakInds = [];
            end
        end

        % Apply threshold value.  Since always finding maxima it will always be
        %   larger than the thresh.
        if ~isempty(thresh)
            m = peakMags>thresh;
            peakInds = peakInds(m);
            peakMags = peakMags(m);
        end



        % Rotate data if needed
        if flipData
            peakMags = peakMags.';
            peakInds = peakInds.';
        end



        % Change sign of data if was finding minima
        if extrema < 0
            peakMags = -peakMags;
            x0 = -x0;
        end
        % Plot if no output desired
        if nargout == 0
            if isempty(peakInds)
                disp('No significant peaks found')
            else
                eval(['figure(''visible'',''' p.displayimage ''');']); 
                plot(1:len0,x0,'.-',peakInds,peakMags,'ro','linewidth',2);
            end
        else
            varargout = {peakInds,peakMags};
        end
    end
    function k3=shortk2(k2)
        k3=[];k3_temp=[];
        if size(k2,2)>2
            for ii=2:size(k2,2)
                if k2(1,ii)-k2(1,ii-1)>3
                    if isempty(k3_temp)==1
                        k3=[k3 k2(:,ii-1)];
                    else
                        k3=[k3 k3_temp];k3_temp=[];
                    end
                else
                    k3_temp=[k3_temp k2(:,ii-1)];
                    [~,k3_tempmax]=max(k3_temp(2,:));k3_temp=k3_temp(:,k3_tempmax);
                    if k2(2,ii)>k3_temp(2,1)
                        k3_temp=k2(:,ii);                    
                    end
                end
            end
            if k2(1,ii)-k2(1,ii-1)>3
                k3=[k3 k2(:,ii)];
            else
                k3=[k3 k3_temp];
            end
        else
            k3=k2;
        end
    end
end
function IM2d_blob=link_neck(IM2d,k2,BWedge,thed2)
% Seperate attached cells into blobs using neck points k2

if isempty(k2)==1
    IM2d_blob=IM2d; % if k2 is empty
elseif size(k2,2)==1
    IM2d_blob=IM2d; % if only one k2
else
    nn=size(k2,2); %number of neck points:nn
    k2points=[k2;zeros(2,nn)]; %if #k2 >1, start connect neck points
    for i=1:nn
        k2points(2:3,i)=[BWedge(k2(1,i),1);BWedge(k2(1,i),2)]; % neck point xy coordinates
    end
    % distance matrix of neck points: k2dist=[starti;endi;dist ..]
    k2dist0=[];
    for i=1:nn-1
        k2dist0=[k2dist0 [i*ones(1,nn-i);i+1:nn]];        
    end     
    k2dist=[k2dist0;pdist(k2points(2:3,:)')];
    % linkage from distance matrix k2dist
    %link_pair=zeros(3,fix(nn/2));
    k2dist=wall_boundary(k2dist,k2,size(BWedge,1)/5);
    link_pair=[];lpend=size(k2dist,2);
    for k=1:lpend
        [a,b]=min(k2dist(3,:));
        if a>thed2
            link_pair=[link_pair k2dist(:,b)];kde=size(link_pair,2);
            for kk=1:2 %delete matrix contains link_pair elements
                k2dist=k2dist(:,k2dist(kk,:)~=link_pair(1,kde));k2dist=k2dist(:,k2dist(kk,:)~=link_pair(2,kde));
            end
        else
            if b==1
                k2dist=k2dist(:,2:end);
            elseif b==size(k2dist,2)
                k2dist=k2dist(:,1:end-1);
            else
                k2dist=[k2dist(:,1:b-1) k2dist(:,b+1:end)];
            end
        end
    end
    IM2d_blob=separate_bob(IM2d,k2points,link_pair);
end
    function IM2d_blob=separate_bob(IM2d,k2points,link_pair)
        for j=1:size(link_pair,2)
            x=[k2points(2,link_pair(1,j)) k2points(2,link_pair(2,j))];
            y=[k2points(3,link_pair(1,j)) k2points(3,link_pair(2,j))];
            if abs(x(2)-x(1))>=abs(y(2)-y(1))
                if x(end)-x(1)>=0
                    X=x(1):x(end);
                else
                    X=x(1):-1:x(end);
                end
                Y=interp1(x,y,X);Y=round(Y);
            else
                if y(end)-y(1)>=0
                    Y=y(1):y(end);
                else
                    Y=y(1):-1:y(end);
                end
                X=interp1(y,x,Y);X=round(X);
            end           
            for jj=1:size(X,2)
                IM2d(Y(1,jj),X(1,jj))=0;
            end
        end
        IM2d_blob=IM2d;
    end
    function k3dist=wall_boundary(k2dist,k2,n)
        k3dist=[];
        for ii=1:size(k2dist,2)
            if abs(k2(1,k2dist(1,ii))-k2(1,k2dist(2,ii)))>n
                k3dist=[k3dist k2dist(:,ii)];
            end
        end
    end
end
function d=distance_fix(M,ma)
n=size(M,1);d=zeros(n,1);
for i=1:n
    d(i,1)=sqrt((M(i,1)-ma(1,1))^2+(M(i,2)-ma(1,2))^2);
end
end

% sub-function
function [IMseries,tinfo]=ImageTiffSeries(tif_name)
outfig=0;

tif_info=imfinfo(tif_name);
Page_num=size(tif_info,1);
for i=1:Page_num
    IMseries{i}=imread(tif_name,i);
    if outfig==1
        eval(['f1=figure(''visible'',''' p.displayimage ''');']); 
        imagesc(IMseries{i});colorbar;colormap jet
        eval([  'print(f1,''-r300'',''-dtiff'',''' 'im_', num2str(i) '.tiff'');'  ]);
    end
end

te=tif_info.Format;tinfo.Format=te;
te=tif_info.Width;tinfo.Width=te;
te=tif_info.Height;tinfo.Height=te;
te=tif_info.BitDepth;tinfo.BitDepth=te;
te=tif_info.ColorType;tinfo.ColorType=te;
te=tif_info.PhotometricInterpretation;tinfo.detail=te;
te=tif_info.XResolution;tinfo.XResolution=te;
te=tif_info.YResolution;tinfo.YResolution=te;
te=tif_info.ResolutionUnit;tinfo.ResolutionUnit=te;
te=tif_info.MaxSampleValue;tinfo.max=te;
te=tif_info.MinSampleValue;tinfo.min=te;

end

% old functions
function Lfix=watershed_fix(L,IM3d)
r=1;
threh=20;

[n,m]=size(L);Lfix=L;
for i=1:n
    for j=1:m
        if L(i,j)==0           
            Lfix(i,j)=water_fix(IM3d,L,i,j,r,threh,n,m);
        end        
    end
end

end
function k2=conv_point2(BWedge,thed,IM2d)
% Obtain neck points k2 from BWedge xy coordinates with thed distance threshold
    n=2;
    [im_a,im_b]=size(IM2d);
    ded_vec=angle_line(BWedge,n);
    k=find(ded_vec<thed);k2=[];
    for iik=1:size(k,1)
        if isempty(find(1:n==BWedge(k(iik),1)))==0||isempty(find(1:n==BWedge(k(iik),2)))==0||isempty(find(im_b-n+1:im_b==BWedge(k(iik),1)))==0||isempty(find(im_a-n+1:im_a==BWedge(k(iik),2)))==0
            k(iik)=0;
        end
    end
    k=k(k~=0);
    if sum(diff(k)./abs(diff(k)))<0 % check order of k equal to BWedge
        k=k(end:-1:1);
    end
    % find neck points on cell boundary: k2
    k0=k;
    %
    for j=1:size(k,1)
       if isempty(find(k-k(j)<5&k-k(j)>0))==0||k(j)~=0
           b=find(k-k(j)<5&k-k(j)>0);
           if isempty(b)==0
               b_list=[[k(j) k(b)'];[ded_vec(k(j)) ded_vec(k(b))']];
               [~,a]=min(b_list(2,:));
               k(k==b_list(1,a))=0;
           end
       end
    end
       %
k2=k0';
    function ded_vec=angle_line(BWedge,n)
    % calculate direction vector of xy in window of boundary
    a=size(BWedge,1);vec=zeros(a,2); %vec: mean vector series on wboundary with window size n
    ded_vec=zeros(a,1); % vecd: dot product between i, i+1 mean vector
    for ii=1:a
        if ii-1<n
            vec(ii,1:2)=mean(diff([BWedge(a-n+ii:a,:);BWedge(1:ii+n,:)]));
        elseif ii>a-n
            vec(ii,1:2)=mean(diff([BWedge(ii-n:a,:);BWedge(1:n-a+ii,:)]));
        else
            vec(ii,1:2)=mean(diff(BWedge(ii-n:ii+n,:)));
        end
        % calculate dplot
        if ii>1
            ded_vec(ii,1)=dot(vec(ii-1,1:2)/sqrt(vec(ii-1,1)^2+vec(ii-1,2)^2),vec(ii,1:2)/sqrt(vec(ii,1)^2+vec(ii,2)^2));
        end
    end
    ded_vec(1,1)=dot(vec(1,1:2)/sqrt(vec(1,1)^2+vec(1,1)^2),vec(a,1:2)/sqrt(vec(a,1)^2+vec(a,1)^2));   
    end
end

