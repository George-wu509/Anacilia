function AnaCilia()
% [Description]



%% == Parameter setting ==
clc;clearvars;
p=parameter_setting();

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
ori_Image = imread(fullFileName,p.rgb); % ori_Image: NxM unit16
warning('off', 'Images:initSize:adjustingMag');
figure;imagesc(ori_Image);

[IMseries,tinfo]=ImageTiffSeries(fullFileName);
[binaryImage,blob]=AnaImage(ori_Image,p);

save('result.mat','IMseries','tinfo','blob');
end
function p=parameter_setting()
p.FileName='WT10%FBS60X1.tif';
p.rgb=4;
p.thresholdValue = 200;

p.se_tophat=30;
p.se_open=15;
p.se_erode=5;
p.num_erode=1;
p.maxgray_o=65536;
p.maxgray_new=4000;
end
function [binaryImage,blob]=AnaImage(IM,p)
% Image process for one image

%% == 1. Image Top-hat filtering ==

originalImage=filter_tophat(IM,p.se_tophat,1);
%figure;imshow(originalImage);


%% == 2. Draw histogram and use threshold value to identify object ==
[~,binaryImage]=imhist_re(originalImage,p.maxgray_o,p.maxgray_new);


%% == 3. Morphologically open image ==
binaryImage=filter_open(binaryImage,p.se_open,1);
%figure;imshow(binaryImage);


%% == 4. Morphologically erode image ==
%binaryImage=filter_erode(binaryImage,p.se_erode,p.num_erode);
%figure;imshow(binaryImage);


%% == 5. Get all the blob properties and draw boundary ==
blob=BlobInfo(binaryImage,originalImage,p);
BlobBounrady(binaryImage,originalImage,p);
eval([ 'save(''result' num2str(p.rgb) '.mat'',''IM'',''binaryImage'',''blob'');'] );

    function [pixelC,binaryImage]=imhist_re(originalImage,maxgray_o,maxgray_new)
        %Draw histogram figure and use threshold value to cut image
        
        [pix, gra] = imhist(originalImage,maxgray_o);
        pixelC=pix(1:maxgray_new,1);grayl=gra(1:maxgray_new,1);
        figure;bar(pixelC);
        hold on;
        maxYValue = ylim;
        line([p.thresholdValue, p.thresholdValue], maxYValue, 'Color', 'r');
        binaryImage = originalImage > p.thresholdValue;
        hold off;
    end
    function [blob,blobMeasurements]=BlobInfo(BW,IM,p)
        % Get all the blob properties and draw boundary
        % BW: binaryImage
        % IM: originalImage
        
        labeledImage = bwlabel(BW, 8);
        blobMeasurements = regionprops(labeledImage, IM, 'all');
        numberOfBlobs = size(blobMeasurements, 1);
        
        figure;imshow(BW); 
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
        
        figure;imshow(IM);
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

% transform
function rgb=image_watershed(BW,BW_obj)
% [ image watershed transform ]
% BW: binary image
% BW_obj=1(object=1), =0(object=0)
% D: distance transform of the complement 
% L: labeling of distance transform of the complement 

    if BW_obj==1
        D=bwdist(~BW);
        D=-D;D(~BW)=-Inf; 
        L=watershed(D);rgb=label2rgb(L);
    elseif BW_obj==0
        D=bwdist(BW);
        D=-D;D(~BW)=-Inf;  
        L=watershed(D);rgb=label2rgb(L);
    end
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

% sub-function
function [IMseries,tinfo]=ImageTiffSeries(tif_name)
outfig=0;

tif_info=imfinfo(tif_name);
Page_num=size(tif_info,1);
for i=1:Page_num
    IMseries{i}=imread(tif_name,i);
    if outfig==1
        f1=figure();set(f1,'Visible','off');
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

