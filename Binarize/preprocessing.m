% =========================================================================
% Vincent Mattana 21128707 2016 MSc IT/Computer Science
% Cosmic Ray Reconstruction
% =========================================================================
clear
%--------------------Paths-------------------------------------------------
path='F:\Ground Truth\Generate\Binarize\';
[filename, pathname] = ...
    uigetfile({'*.bmp;*.jpg;',...
    'Picture Files (*.bmp,*.jpg)';},'Select source image for processing');

%Input Files, and details thereof
%--------------------------------------------------------------------------
input=imread(strcat(strcat(pathname,'\'),filename) );
inputGREY=rgb2gray(double(input));
inputHSV=rgb2hsv(double(input));
inputBW = im2bw(input, 0.75);

dimension=size(inputHSV);
Rdim=dimension(1);
Cdim=dimension(2);
thresh=0.5;
HOUR_MARKER_COLOUR =1;

input_message=['Input image has ',num2str(Rdim),' rows, and ',num2str(Cdim),' columns.'];
disp(input_message)

% % Construct a questdlg with three options
% choice = questdlg('What colour are the hour markers?', ...
% 	'C O S M I C U S 2 0 1 6', ...
% 	'W H I T E','B L A C K','X');
% % Handle response
% switch choice
%     case 'White'
%         disp([choice ' hour markers.'])
%         HOUR_MARKER_COLOUR_temp = 1;
%     case 'Black'
%         disp([choice ' hour markers.'])
%         HOUR_MARKER_COLOUR_temp = 0;
% end

%Sprockets - Finding via Correlation
%--------------------------------------------------------------------------
idealSprockets = sprockets(input);

CC_Sprockets = bwconncomp(idealSprockets);
sprocket_message=[num2str(CC_Sprockets.NumObjects),' sprockets found.'];
disp(sprocket_message)

%Sprocket Removal
%--------------------------------------------------------------------------
desprocketed_TEST=inputHSV;
for y=1:Cdim,
    for x=1:Rdim,
        if (idealSprockets(x,y) == 1)
            desprocketed_TEST(x,y,:) = NaN;
        end
    end
end
desprocketed=desprocketed_TEST(:,:,3)/255;
%imwrite(desprocketed,strcat(path,'01-desprocketed.bmp'),'bmp');

%Segmentation
%--------------------------------------------------------------------------

binarized_0=zeros(Rdim,1);
idealHM_0=zeros(Rdim,1);

for i=0:18
    %---------------Cropping Window Parameters---------------------------
    %Binary =blockproc( reshape(1:18,315,Rdim)', [1,1], @(x) BI(:,:,x.data), 'UseParallel',true);
    
    
    window = [i*(round(Cdim/19)),0,(round(Cdim/19)-1),Rdim];%[i*316,0,315,Rdim]
    HSV_segment = imcrop(desprocketed, window);
    ideal_SL_segment=false(size(HSV_segment));
    ideal_HM_segment=false(size(HSV_segment));
    
    dimension_segment=size(HSV_segment);
    Rdim_seg_0=dimension_segment(1);
    Cdim_seg_0=dimension_segment(2);
    
    %---------------Thresholding Algorithms---------------------------
    binarized_mini_0=zeros(Rdim,(round(Rdim/5)));
    fifths = [0,round(Rdim/5),round(2*(Rdim/5)),round(3*(Rdim/5)),round(4*(Rdim/5)),Rdim];
    
    for j=1:5
        %---------------Cropping Window Parameters---------------------------
        window_2 = [0,fifths(j),Cdim,(fifths(j+1)-fifths(j))-1];
        mini_segment = imcrop(HSV_segment, window_2);
        
        
        mini_dimension_segment=size(mini_segment);
        ydim_mini_seg=mini_dimension_segment(1);
        xdim_mini_seg=mini_dimension_segment(2);
        
        thresh = 0.5;
        binarized_mini_Segment = im2bw(mini_segment,thresh);
        for k=1:250
            whitecount = sum(sum(binarized_mini_Segment == 1))/(ydim_mini_seg*xdim_mini_seg);
            if whitecount > 0.22
                thresh = thresh+0.002;
            end
            if whitecount < 0.18
                thresh = thresh-0.001;
            end
            binarized_mini_Segment = im2bw(mini_segment, thresh);
        end
        
        %---------------Concat Segments---------------------------
        if(j==1)
            binarized_mini_0=binarized_mini_Segment;
        end
        if(j>1)
            binarized_mini_0=cat(1,binarized_mini_0,binarized_mini_Segment);
        end
    end
    %figure, imshow(binarized_mini_0)
    binarized_section=binarized_mini_0;
    
    %---------------Removing Noise---------------------------
    binarized_section = imerode(binarized_section, strel('rectangle', [2 2]));
    binarized_section = imdilate(binarized_section, strel('rectangle', [2 2]));
    binarized_section = imerode(binarized_section, strel('rectangle', [2 2]));
    binarized_section = imdilate(binarized_section, strel('rectangle', [2 2]));
    binarized_section = imerode(binarized_section, strel('rectangle', [2 2]));
    
    %---------------Getting Hour Markers--------------------------
    
    if(HOUR_MARKER_COLOUR==1)
        totalsHM=sum(binarized_section,1,'omitnan')/dimension_segment(1);
        for x=1:dimension_segment(2)-1,
            for y=1:dimension_segment(1),
                if (totalsHM(x) > (0.75 *max(totalsHM)))
                    ideal_HM_segment(y, x) = 1;
                end
            end
        end
    end
    if(HOUR_MARKER_COLOUR==0)
        HSV_segment_comp=imcomplement(HSV_segment);
        binarized_section_comp=im2bw(HSV_segment_comp,0.5);
        totalsHM=sum(binarized_section_comp,1,'omitnan')/dimension_segment(1);
        for x=1:dimension_segment(2)-1,
            for y=1:dimension_segment(1),
                if (totalsHM(x) > (0.75 *max(totalsHM)))
                    ideal_HM_segment(y, x) = 1;
                end
            end
        end
    end
    
    
    %---------------Dilating hour markers--------------------------
    hmWidth=sum(ideal_HM_segment,2,'omitnan')/dimension_segment(2);
    if (hmWidth<16)
        ideal_HM_segment=imdilate(ideal_HM_segment,strel('rectangle', [1 16]));
    end
    
    %     data_segment=binarized_section;
    %     for y=1:dimension_segment(2)-1,
    %         for x=1:dimension_segment(1)-1,
    %             if (ideal_HM_segment(x,y) == 1)
    %             data_segment(x,y)=0;
    %             end
    %             if (ideal_SL_segment(x,y) == 1)
    %             data_segment(x,y)=0;
    %             end
    %         end
    %     end
    %
    %     ideal_data=false(size(data_segment));
    %     data_segment=imclose(data_segment,strel('diamond',3));
    %     CC = bwconncomp(data_segment);
    %     numPixels = cellfun(@numel,CC.PixelIdxList);
    %     [biggest,idx] = max(numPixels);
    %     ideal_data(CC.PixelIdxList{idx}) = 1;
    %     %CC.PixelIdxList(idx)=[];
    %
    %     data_segment=ideal_data;
    
    %     %---------------Finding hour markers via Hough--------------------------
    %     start_angle = -2.5; % >-90 && <90
    %     end_angle = 2.5; % >-90 && <90
    %     theta_resolution = 0.75;
    %     [H T R] = hough(binarized_section,'RhoResolution',0.9, 'Theta', start_angle:theta_resolution:end_angle);
    %     P  = houghpeaks(H,1,'threshold',ceil(max(H(:))));
    %     lines = houghlines(binarized_section, T, R, P, 'MinLength',3);
    %
    %     deHourMarked_section = descaled_section;
    %     for k=1:length(lines);
    %         %---------------Removing Hour Markers via Bresenham Algorithm-------------
    %         [x, y] = bresenham(lines(k).point1(1), lines(k).point1(2), ...
    %                            lines(k).point2(1), lines(k).point2(2));
    %         hm(y, x) = 1;
    %     end
    %
    %     for x=1:dimension_segment(2)-1,
    %         for y=1:dimension_segment(1),
    %             if(hm(y, x)==1)
    %                 deHourMarked_section(y, x) = 0;
    %             end
    %         end
    %     end
    %     %---------------Display Found Lines---------------------------
    %     figure, imshow(descaled_section), hold on
    %     for k = 1:length(lines)
    %         xy = [lines(k).point1; lines(k).point2];
    %         plot(xy(:,1), xy(:,2), 'g.-', 'LineWidth',2);
    %     end
    %     hold off
    
    %---------------Concat Segments---------------------------
    if(i==0)
        binarized_0=binarized_section;
        idealHM_0=ideal_HM_segment;
    else
        binarized_0=cat(2,binarized_0,binarized_section);
        idealHM_0=cat(2,idealHM_0,ideal_HM_segment);
    end
    segment_message=['Segment_',num2str(i),' processed.'];
    disp(segment_message)
end
imwrite(binarized_0,strcat(path,'02-Binarized.bmp'),'bmp');
imwrite(idealHM_0,strcat(path,'05c-idealHM.bmp'),'bmp');

CC_HM = bwconncomp(idealHM_0);
hm_message=[num2str(CC_HM.NumObjects),' hour markers found.'];
disp(hm_message)

%ideal compilation
%--------------------------------------------------------------------------
dimension_ideal=size(idealHM_0);
Rdim_ideal=dimension_ideal(1);
Cdim_ideal=dimension_ideal(2);
ideal=false(size(inputBW));
for y=1:Cdim_ideal,
    for x=1:Rdim_ideal,
        if ((idealSprockets(x,y) == 1)||(idealHM_0(x,y)==1))
            ideal(x,y)=1;
        end
    end
end
imwrite(ideal,strcat(path,'06-GroundTruth.bmp'),'bmp');

overlay=input;
for y=1:Cdim_ideal,
    for x=1:Rdim_ideal,
        if (idealSprockets(x,y) == 1)
            overlay(x,y,:)=NaN;
        end
        %         if (idealSL_0(x,y)==1)
        %         overlay(x,y,1)=NaN;
        %         end
        if (idealHM_0(x,y)==1)
            overlay(x,y,3)=255;
        end
    end
end
imwrite(overlay,strcat(path,'01-Overlay.bmp'),'bmp');
disp('First pass complete.')

%PASS TWO (2)
%--------------------------------------------------------------------------

desprocketed_v2=desprocketed;

%---------------counting hm---------------------------
count=0;
totalCount=0;
hmMiddleEND=(CC_HM.NumObjects);
hmMiddle = zeros(1,hmMiddleEND);
for y=2:Cdim_ideal,
    if(idealHM_0(1,y)==1)
        count=count+1;
    end
    if(idealHM_0(1,y)==0)&&(idealHM_0(1,y-1)==1)
        totalCount=totalCount+1;
        hmMiddle(1,totalCount)=y-round(count/2);
        count=0;
    end
end

hmMiddle=[hmMiddle,Cdim_ideal];

for i=1:size(hmMiddle,2)-1
    %---------------Cropping Window Parameters 2---------------------------
    window = [hmMiddle(i),0,(hmMiddle(i+1)-hmMiddle(i)),Rdim];
    idealHM_0_Segment = imcrop(idealHM_0, window);
    input_v2_Segment=imcrop(input, window);
    binarized_Segment=imcrop(binarized_0, window);
    desprocketed_v2_Segment = imcrop(desprocketed_v2, window);
    idealSprockets_v2_Segment = imcrop(idealSprockets, window);
    
    
    dimension_segment=size(binarized_Segment);
    Rdim_seg=dimension_segment(1);
    Cdim_seg=dimension_segment(2);
    
    %     %---------------Removing found ideal lines-------------
    %     for x=1:Cdim_seg-1,
    %         for y=1:Rdim_seg,
    %             if(ideal_SL_segment(y, x)==1)
    %                 descaled_section(y, x) = 0;
    %             end
    %         end
    %     end
    %
    %     %---------------Tiding Up after removal---------------------------
    %     descaled_section = imerode(descaled_section, strel('rectangle', [4 1]));
    %     descaled_section = imerode(descaled_section, strel('rectangle', [1 4]));
    
    
    
    %   %---------------Scale lines 2---------------------------
    %     mask=imcomplement(idealSL_0_Segment);
    %     sl_Segment=desprocketed_v2_Segment;
    %     sl_Segment_v2=desprocketed_v2_Segment;
    %     for x=1:Rdim_seg-1,
    %         for y=1:Cdim_seg-1,
    %             if (mask(x,y) == 1)
    %             sl_Segment(x,y,:)=NaN;
    %             end
    %         end
    %     end
    %     sl_Segment=im2bw(sl_Segment,thresh-(thresh*0.2));
    %     sl_Segment=bwmorph(sl_Segment,'skel',Inf);
    %     sl_Segment=bwmorph(sl_Segment, 'spur', 10);
    %     sl_Segment=imfill(sl_Segment,8,'holes');
    %
    %     sl_Segment=imdilate(sl_Segment,strel('rectangle', [2 16]));
    %     sl_Segment=imfill(sl_Segment,8,'holes');
    %     sl_Segment=bwmorph(sl_Segment,'skel',Inf);
    %     sl_Segment=bwmorph(sl_Segment, 'spur', 10);
    %
    %     sl_Segment=imdilate(sl_Segment,strel('rectangle', [2 8]));
    %     sl_Segment=imfill(sl_Segment,8,'holes');
    %
    %     sl_Segment=imdilate(sl_Segment,strel('rectangle', [2 4]));
    %
    %     sl_Segment=bwmorph(sl_Segment,'skel',Inf);
    %     sl_Segment=bwmorph(sl_Segment, 'spur', 20);
    %
    %     sl_Segment=imdilate(sl_Segment,strel('octagon', 3));
    
    %---------------Scale lines 2---------------------------
    ideal_SL_segment=scalelines(binarized_Segment);
    ideal_SL_segment=imdilate(ideal_SL_segment,strel('octagon', 3));
    ideal_SL_segment=bwmorph(ideal_SL_segment,'skel',Inf);
    ideal_SL_segment=bwmorph(ideal_SL_segment, 'spur', 20);
    ideal_SL_segment=imdilate(ideal_SL_segment,strel('octagon', 3));
    
    %---------------hour markers 2---------------------------
    mask=imcomplement(idealHM_0_Segment);
    hm_Segment=desprocketed_v2_Segment;
    for x=1:Rdim_seg,
        for y=1:Cdim_seg-1,
            if (mask(x,y) == 1)
                hm_Segment(x,y,:)=NaN;
            end
        end
    end
    hm_Segment=im2bw(hm_Segment,thresh);
    hm_Segment=bwmorph(hm_Segment, 'spur', 25);
    hm_Segment=imfill(hm_Segment,8,'holes');
    hm_Segment = bwareaopen(hm_Segment, 25);
    %---------------Data Line Approx---------------------------
    data_Segment=binarized_Segment;
    for y=1:Cdim_seg,
        for x=1:Rdim_seg,
            if ((ideal_SL_segment(x,y)==1)||(hm_Segment(x,y)==1)||(idealSprockets_v2_Segment(x,y)==1))
                data_Segment(x,y)=0;
            end
        end
    end
    
    data_Segment= bwareaopen(data_Segment,25);
    %     data_Segment=edge(data_Segment,'canny');
    %     data_Segment= bwareaopen(data_Segment,30);
    data_Segment=imdilate(data_Segment,strel('octagon', 6));
    data_Segment=bwmorph(data_Segment,'skel',Inf);
    data_Segment=bwmorph(data_Segment, 'spur', 200);
    %     data_Segment= bwareaopen(data_Segment,30);
    data_Segment=imdilate(data_Segment,strel('octagon', 6));
    CC = bwconncomp(data_Segment);
    
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    data_Segment=zeros(size(data_Segment));
    data_Segment(CC.PixelIdxList{idx}) = 1;
    figure, imshow(data_Segment)
    
    %---------------Concat Segments---------------------------
    if(i==1)
        idealSL_v2=ideal_SL_segment;
        idealHM_v2=hm_Segment;
        input_v2=input_v2_Segment;
        idealSprockets_v2=idealSprockets_v2_Segment;
        data=data_Segment;
    else
        idealSL_v2=cat(2,idealSL_v2,ideal_SL_segment);
        idealHM_v2=cat(2,idealHM_v2,hm_Segment);
        input_v2=cat(2,input_v2,input_v2_Segment);
        idealSprockets_v2=cat(2,idealSprockets_v2,idealSprockets_v2_Segment);
        data=cat(2,data,data_Segment);
    end
    
    segment_message=['Segment_',num2str(i),' processed, again.'];
    disp(segment_message)
end

dimension_v2=size(idealSL_v2);
Rdim_v2=dimension_v2(1);
Cdim_v2=dimension_v2(2);

%ideal compilation
%--------------------------------------------------------------------------
ideal_v2=zeros(size(idealSL_v2));
for y=1:Cdim_v2,
    for x=1:Rdim_v2,
        if ((idealSL_v2(x,y)==1)||(idealHM_v2(x,y)==1)||(idealSprockets_v2(x,y)==1))
            ideal_v2(x,y)=1;
        end
    end
end
imwrite(ideal_v2,strcat(path,'ideal_2.bmp'),'bmp');

overlay_v2=input_v2;
for y=1:Cdim_v2-1,
    for x=1:Rdim_v2-1,
        if (idealSprockets_v2(x,y) == 1)
            overlay_v2(x,y,:)=NaN;
        end
        if (idealSL_v2(x,y)==1)
            overlay_v2(x,y,:)=NaN;
        end
        if (idealHM_v2(x,y)==1)
            overlay_v2(x,y,:)=1;
        end
        if (data(x,y)==1)
            overlay_v2(x,y,:)=1;
        end
    end
end
imwrite(overlay_v2,strcat(path,'02-Overlay.bmp'),'bmp');
disp('Second pass complete.')

%figure, imshow(inputGREY)


%**************************************************************************
disp('Operation complete.')
%**************************************************************************