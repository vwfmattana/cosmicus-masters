%   ???????????  Vincent Mattana  -??????????
%    ,.-~*����`*�~-.�-(Cosmicus)-,.-~*����`*�~-.�
%     ,.-~*����`*�~-.�-(M.Sc)-,.-~*����`*�~-.�
%     ,.-~*����`*�~-.�-(2016)-,.-~*����`*�~-.�
clear;
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
disp('       _        _        _        _        _        _           _        _        _        _        _        _')
disp('      / /\     / /\     / /\     / /\     / /\     / /\       /\ \     /\ \     /\ \     /\ \     /\ \     /\ \')
disp('     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \     /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \')
disp('    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \   / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \')
disp('   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/   \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \')
disp('  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \  Welcome  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \')
disp('  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / /  _     _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /')
disp('   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\   /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /')
disp('    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /   \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /')
disp('     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /     \  / /   \  / /   \  / /   \  / /   \  / /   \  / /')
disp('      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/       \/_/     \/_/     \/_/     \/_/     \/_/     \/_/')
%'
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
disp('Operation Initiated.')
path='/home/vincent/cosmicus/code/Binarize';
[filename, pathname] = ...
    uigetfile({'*.bmp;*.jpg;',...
    'Picture Files (*.bmp,*.jpg)';},'Select ORIGINAL image for processing');
input=imread(strcat(strcat(pathname,'/'),filename) );

path='/home/vincent/cosmicus/code/Binarize';
[filename, pathname] = ...
    uigetfile({'*.bmp;*.jpg;',...
    'Picture Files (*.bmp,*.jpg)';},'Select USER image for processing');
inputUSER=imread(strcat(strcat(pathname,'/'),filename) );
%   ??????????? Input Files, and details thereof -??????????

inputGREY=rgb2gray(double(input));
inputHSV=rgb2hsv(double(input));
inputBW = im2bw(input, 0.5);
dimension=size(inputGREY);
Rdim=dimension(1);
Cdim=dimension(2);

% % % Shortcut
% path='F:\Ground Truth\Generate\Binarize\';
% [filename, pathname] = ...
%     uigetfile({'*.bmp;*.jpg;',...
%     'Picture Files (*.bmp,*.jpg)';},'Select binarized image:');
% binary=imread(strcat(strcat(pathname,'\'),filename) );

HOUR_MARKER_COLOUR=0;
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
disp('       _        _        _        _        _        _             _        _        _        _        _        _')
disp('      / /\     / /\     / /\     / /\     / /\     / /\         /\ \     /\ \     /\ \     /\ \     /\ \     /\ \')
disp('     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \       /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \')
disp('    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \     / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \')
disp('   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/     \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \')
disp('  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \  Sprocket   / / /     \ \ \   / / /     \ \ \   / / /     \ \ \')
disp('  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / / _Removal _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /')
disp('   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\     /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /')
disp('    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /     \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /')
disp('     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /       \  / /   \  / /   \  / /   \  / /   \  / /   \  / /')
disp('      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/         \/_/     \/_/     \/_/     \/_/     \/_/     \/_/')
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%histeq(input);
idealSprockets = sprockets(input);
imwrite(idealSprockets,strcat(path,'idealSprockets.bmp'),'bmp');
desprocketed_input=input;
desprocketed_inputUSER=inputUSER;
for y=1:Cdim
    for x=1:Rdim
        if (idealSprockets(x,y) == 1)
            desprocketed_input(x,y,:) = NaN;
            desprocketed_inputUSER(x,y,:) = NaN;
        end
    end
end
%desprocketing_mask_HSV=rgb2hsv(double(desprocketed_input));
%desprocketed=desprocketing_mask_HSV(:,:,3)/255;
disp('Sprockets Removed.')

%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
disp('       _        _        _        _        _        _             _        _        _        _        _        _')
disp('      / /\     / /\     / /\     / /\     / /\     / /\         /\ \     /\ \     /\ \     /\ \     /\ \     /\ \')
disp('     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \       /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \')
disp('    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \     / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \')
disp('   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/     \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \')
disp('  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \  Binar-     / / /     \ \ \   / / /     \ \ \   / / /     \ \ \')
disp('  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / /  ization _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /')
disp('   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\     /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /')
disp('    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /     \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /')
disp('     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /       \  / /   \  / /   \  / /   \  / /   \  / /   \  / /')
disp('      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/         \/_/     \/_/     \/_/     \/_/     \/_/     \/_/')
%'
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�

fun0=@(block_struct) rgb2gray(block_struct.data); 'BorderSize'; [5 5];
gray =blockproc( desprocketed_input, [10,10], fun0);
grayUSER =blockproc( desprocketed_inputUSER, [10,10], fun0);
%gray2=histeq(gray);
%imwrite(gray,strcat(path,'01-Desprocketed(Grayscale-eq).bmp'),'bmp');

%Long Road
fun1=@(block_struct) binarize(block_struct.data); 'BorderSize'; [5 5];
[mb,nb] = bestblk([Rdim Cdim],50);
binary =blockproc( grayUSER, [mb,nb], fun1);
imwrite(binary,strcat(path,'02-Binary.bmp'),'bmp');

%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
disp('       _        _        _        _        _        _             _        _        _        _        _        _')
disp('      / /\     / /\     / /\     / /\     / /\     / /\         /\ \     /\ \     /\ \     /\ \     /\ \     /\ \')
disp('     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \       /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \')
disp('    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \     / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \')
disp('   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/     \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \')
disp('  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \  Hour -     / / /     \ \ \   / / /     \ \ \   / / /     \ \ \')
disp('  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / /  Markers _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /')
disp('   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\     /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /')
disp('    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /     \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /')
disp('     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /       \  / /   \  / /   \  / /   \  / /   \  / /   \  / /')
disp('      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/         \/_/     \/_/     \/_/     \/_/     \/_/     \/_/')
%'
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�


fun2=@(block_struct) hourmarkers(block_struct.data,HOUR_MARKER_COLOUR); 'BorderSize';
hm =blockproc( grayUSER, [Rdim,floor(Cdim/20)], fun2);
%imwrite(hm,strcat(path,'05-hm.bmp'),'bmp');

totalHMCount=0;
current_gap=[];
gap_count=0;
width_count=0;
for y=2:size(hm,2)-1
    if(hm(1,y)==1) %White Col.
        if(hm(1,y-1)==0) %Prev is white
            current_gap=[current_gap, gap_count];
            gap_count=0;
        end
        width_count=width_count+1;
    end
    
    if(hm(1,y)==0)  %Black Col
        if(hm(1,y-1)==1) %Prev is white
            totalHMCount=totalHMCount+1;
            width_count=0;
        end
        gap_count=gap_count+1;
    end
end
mean_hm_gapSize=round(mean2(current_gap));
median_hm_gapSize=round(median(current_gap));


gap_count=1;
width_count=0;
totalHMCount=0;
FIRST_HM_FOUND=0;
hmMiddleEND=35;
hmMiddle = zeros(1,hmMiddleEND-1);
gap_array=[];
for y=2:Cdim
    if(FIRST_HM_FOUND==0)%first hm
        if(hm(1,y)==0) && (hm(1,y-1)==1)   %Black Col, Prev is white
            totalHMCount=totalHMCount+1;
            hmMiddle(1,totalHMCount)=y;%-round(width_count/2);
            width_count=0;
            FIRST_HM_FOUND=1;
        end
        if(hm(1,y)==1) %White Col.
            width_count=width_count+1;
        end
    else
        if(hm(1,y)==1) %White Col.
            if((gap_count<mean_hm_gapSize-20)&&(gap_count>0))%Gap is present & smaller than median gap size
                hm(:,y)=0;%delete col.
            else
                gap_count=0;
                width_count=width_count+1;
            end
        end
        
        if(hm(1,y)==0)  %Black Col
            if(hm(1,y-1)==1) %Prev is white
                totalHMCount=totalHMCount+1;
                hmMiddle(1,totalHMCount)=y-round(width_count/2);
                width_count=0;
                gap_count=0;
            end
            gap_count=gap_count+1;
        end
    end
    
end

hmMiddle=[hmMiddle,Cdim];
toRemove=[];
for i=1:size(hmMiddle,2)
    if hmMiddle(1,i) ==0
        toRemove = [toRemove; i];
    end
end

hmMiddle=hmMiddle(setdiff(1:length(hmMiddle),toRemove));
hmMiddle=[0,hmMiddle];
ideal_HM=false(size(hm));
% ---------------------Tilting-------------------------
for i=2:size(hmMiddle,2)-1
    % -----------------Making Window----------------------
    window_HM = [hmMiddle(i)-40,0,80,Rdim];
    hm_gray_line_Segment = imcrop(grayUSER, window_HM);
    hm_line_Segment=imcrop(hm, window_HM);
    HM_Segment_comp=imcomplement(hm_gray_line_Segment);
    
    %HM_Segment_comp=im2bw(HM_Segment_comp,0.5);
    HM_edge = edge(hm_gray_line_Segment,'canny');
    
    %   ??????????? Finding hour markers via Hough--------------------------
    start_angle = -7; % >-90 && <90
    end_angle = 7; % >-90 && <90
    theta_resolution = 0.5;
    HM_HOUGH_THRESH=0.75;
    [H, T, R] = hough(HM_edge,'RhoResolution',2, 'Theta', start_angle:theta_resolution:end_angle);
    P  = houghpeaks(H,1,'threshold',ceil(HM_HOUGH_THRESH*max(H(:))));
    lines_hm = houghlines(HM_edge, T, R, P);
    
    %# shearing transforma
    slopes = vertcat(lines_hm.point2) - vertcat(lines_hm.point1);
    slopes = slopes(:,2) ./ slopes(:,1);

    w = atand(Inf) - atand(slopes(1));
    delta = round(Rdim*tan(deg2rad(w))/2);
    [x, y] = bresenham(hmMiddle(i)-delta, 1, hmMiddle(i)+delta, Rdim);
    ind = sub2ind( size(ideal_HM), y, x );
    ideal_HM(ind) = 1;
    hm=imdilate(ideal_HM,strel('octagon', 6));
end
ideal_HM_shifted=false(size(hm));
for i=2:size(hmMiddle,2)-1
    candidate_col=[];
    candidate_totals=[];
    for k=-15:15
        % -----------------Making Window----------------------
        window_HM_still = [hmMiddle(i)-15,0,31,Rdim];
        window_HM = [hmMiddle(i)+k,0,31,Rdim];
        hm_line_Segment=imcrop(hm, window_HM_still);
        gray_line_Segment=imcrop(grayUSER, window_HM);
        %binary_line_Segment=imcrop(gray, window_HM);
        
        
        % ---------------------Cropping Grayscale with CC-mask -------------------------
        hm_only=gray_line_Segment;
        for y=1:size(hm_line_Segment,2)
            for z=1:Rdim
                if(hm_line_Segment(z,y)==0)
                    hm_only(z,y)=NaN;
                end
                
            end
        end
        totalsHM=sum(sum(hm_only,'omitnan'));
        candidate_col=[candidate_col hmMiddle(i)+k];%binarized_Segment
        candidate_totals=[candidate_totals totalsHM];
        
        % :}  <- Skillie
        
    end
    [M,INDEX] = min(candidate_totals);
    loc_hm=candidate_col(INDEX);
    
    ideal_HM_shifted(1:Rdim,loc_hm:loc_hm+31)=hm_line_Segment;
end
%jiggling HM
hm=ideal_HM_shifted;

gap_count=1;
width_count=0;
totalHMCount=0;
FIRST_HM_FOUND=0;
hmMiddleEND=35;
hmMiddle = zeros(1,hmMiddleEND-1);
gap_array=[];
for y=2:Cdim
    if(FIRST_HM_FOUND==0)%first hm
        if(hm(1,y)==1) && (hm(1,y-1)==0)   %Black Col, Prev is white
            totalHMCount=totalHMCount+1;
            hmMiddle(1,totalHMCount)=y;
        end
    end
end

hmMiddle=[hmMiddle,(Cdim+1)];
toRemove=[];
for i=1:size(hmMiddle,2)
    if hmMiddle(1,i) ==0
        toRemove = [toRemove; i];
    end
end

hmMiddle=hmMiddle(setdiff(1:length(hmMiddle),toRemove));
hmMiddle=[0,hmMiddle];

%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
disp('       _        _        _        _        _        _             _        _        _        _        _        _')
disp('      / /\     / /\     / /\     / /\     / /\     / /\         /\ \     /\ \     /\ \     /\ \     /\ \     /\ \')
disp('     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \       /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \')
disp('    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \     / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \')
disp('   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/     \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \')
disp('  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \    Segment  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \')
disp('  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / / Creation _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /')
disp('   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\     /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /')
disp('    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /     \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /')
disp('     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /       \  / /   \  / /   \  / /   \  / /   \  / /   \  / /')
disp('      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/         \/_/     \/_/     \/_/     \/_/     \/_/     \/_/')
%'
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
data_values = [];
deltaN_values =[];
avg_data_Values = [];
percent_Increase_Values=[];
for i=1:size(hmMiddle,2)-1
    % ---------------------Making Window-------------------------
    window = [hmMiddle(i),0,(hmMiddle(i+1)-hmMiddle(i))-1,Rdim];
    hm_Segment = imcrop(hm, window);
    binarized_Segment=imcrop(binary, window);
    
    gray_Segment=imcrop(gray, window);
    gray_SegmentUSER=imcrop(grayUSER, window);
    input_Segment=imcrop(input, window);
    desprocketed_Segment = imcrop(desprocketed_input, window);
    %     idealSprockets_v2_Segment = imcrop(idealSprockets, window);
    
    dimension_segment=size(binarized_Segment);
    Rdim_seg=dimension_segment(1);
    Cdim_seg=dimension_segment(2);
    
    % ---------------------Scalelines.m Blockproc-------------------------
    fun3=@(block_struct) scalelines(block_struct.data); 'BorderSize';
    sl =blockproc( binarized_Segment, [round(Rdim_seg/60),Cdim_seg], fun3);
    
    dimension_segment_sl=size(sl);
    Rdim_seg_SL=dimension_segment_sl(1);
    Cdim_seg_SL=dimension_segment_sl(2);
    
    %=================================== Jiggling==========================
    height_count=0;
    candidate_rows=[];
    candidate_locations=[];
    lenient_mean_gapSize=12;
    sl_temp=zeros(size(sl));
    for y=2:Rdim_seg_SL
        
        if(sl(y,1)==1) %White Col.
            height_count=height_count+1;
        end
        if(sl(y,1)==0)  %Black Col
            if(sl(y-1,1)==1)%Prev is white
                if(y>12) && (y<Rdim_seg_SL-12) && (y>round(height_count/2))%not near borders
                    for j=-9:9
                        candidate_rows=[candidate_rows sum(binarized_Segment((y-1-round(height_count/2))+j,:))];%binarized_Segment
                        candidate_locations=[candidate_locations (y-1-round(height_count/2))+j];
                    end
                    
                    for r=1:9
                        toRemove=[];
                        for s=1:size(candidate_rows,2)-1
                            if ((candidate_locations(1,s+1)-candidate_locations(1,s)<lenient_mean_gapSize))
                                if candidate_rows(1,s+1) > candidate_rows(1,s)
                                    toRemove = [toRemove; s];
                                else
                                    toRemove = [toRemove; s+1];
                                end
                            end
                        end
                        candidate_rows=candidate_rows(setdiff(1:length(candidate_rows),toRemove));
                        candidate_locations=candidate_locations(setdiff(1:length(candidate_locations),toRemove));
                    end
                    
                    [~,INDEX] = max(candidate_rows);
                    loc=candidate_locations(1,INDEX);
                    %sl(y-1:y-1-height_count,:)=0;
                    
                    toRemove=[];
                    for s=-7:7
                        if (loc>7)&&(loc<Rdim_seg_SL-7) && (sl_temp(loc+s,1)==1)
                            toRemove = [toRemove; INDEX];
                            break;
                        end
                    end
                    candidate_rows=candidate_rows(setdiff(1:length(candidate_rows),toRemove));
                    candidate_locations=candidate_locations(setdiff(1:length(candidate_locations),toRemove));
                    
                    [M,INDEX] = max(candidate_rows);
                    loc=candidate_locations(1,INDEX);
                    sl_temp(loc:loc,:)=1;
                    height_count=0;
                    candidate_rows=[];
                    candidate_locations=[];
                end
            end
        end
        
    end
    sl=sl_temp;
    %sl_temp=imdilate(sl_temp,strel('rectangle', [3 3]));
    % ---------------------Getting average gap size-------------------------
    gap_count=1;
    height_count=0;
    totalSLCount=0;
    current_gap=[];
    SL_Locations = zeros(1,65);
    for y=2:Rdim_seg_SL
        if(sl(y,1)==1) %White Col.
            if(sl(y-1,1)==0) %Prev is white
                current_gap=[current_gap, gap_count];
                gap_count=0;
            end
            height_count=height_count+1;
        end
        
        if(sl(y,1)==0)  %Black Col
            if(sl(y-1,1)==1) %Prev is white
                totalSLCount=totalSLCount+1;
                SL_Locations(1,totalSLCount)=y-round(height_count/2);
                height_count=0;
            end
            gap_count=gap_count+1;
        end
    end
    mean_gapSize=median(current_gap);
    gap_count=0;
    % ---------------------Placing missing scale lines, based on avg x10-------------------------
    for u=1:10
        for y=2:Rdim_seg_SL
            if(sl(y,1)==1) %White Col.
                if(sl(y-1,1)==0) %Prev is black
                    if(gap_count > mean_gapSize)&&(y>gap_count)
                        start=round(y-gap_count+mean_gapSize-1);
                        fin=round(y-gap_count+mean_gapSize);
                        sl(start:fin,:)=1;
                    end
                    gap_count=0;
                end
                
            end
            totalsSL=sum(sl,2,'omitnan');
            moving_total=0;
            if(sl(y,1)==0) %Black Col
                gap_count=gap_count+1;
                totalsSL=sum(gray_Segment,2,'omitnan');
            end
        end
        gap_count=0;
    end
    % ---------------------Getting average gap size-------------------------
    gap_count=1;
    height_count=0;
    totalSLCount=0;
    current_gap=[];
    SL_Locations = zeros(1,65);
    for y=2:Rdim_seg_SL
        if(sl(y,1)==1) %White Col.
            if(sl(y-1,1)==0) %Prev is white
                current_gap=[current_gap, gap_count];
                gap_count=0;
            end
            height_count=height_count+1;
        end
        
        if(sl(y,1)==0)  %Black Col
            if(sl(y-1,1)==1) %Prev is white
                totalSLCount=totalSLCount+1;
                SL_Locations(1,totalSLCount)=y-round(height_count/2);
                height_count=0;
            end
            gap_count=gap_count+1;
        end
    end
    mean_gapSize=median(current_gap);
    % ---------------------Deleting needless SL-------------------------
    gap_count=1;
    height_count=0;
    totalSLCount=0;
    current_gap=[];
    average_gap=[];
    FIRST_SL_FOUND=0;
    SL_Locations = zeros(1,65);
    for y=2:Rdim_seg_SL
        if(FIRST_SL_FOUND==0)%first sl
            if(sl(y,1)==0) && (sl(y-1,1)==1)   %Black Col
                totalSLCount=totalSLCount+1;
                SL_Locations(totalSLCount,1)=y;%-round(width_count/2);
                height_count=0;
                FIRST_SL_FOUND=1;
            end
            if(sl(y,1)==1) %White Col.
                height_count=height_count+1;
            end
        else
            if(sl(y,1)==1) %White Col.
                if(gap_count<13)&&(gap_count>0)%Gap is present
                    sl(y,:)=0;
                else
                    gap_count=0;
                    height_count=height_count+1;
                end
                
            end
        end
        if(sl(y,1)==0)  %Black Col
            if(sl(y-1,1)==1) %Prev is white
                totalSLCount=totalSLCount+1;
                SL_Locations(1,totalSLCount)=y;%-round(width_count/2);
                height_count=0;
            end
            gap_count=gap_count+1;
        end
    end
    sl=imdilate(sl,strel('rectangle', [3 3]));
    gap_count=0;
    %=================================== Jiggling==========================
    height_count=0;
    candidate_rows=[];
    candidate_locations=[];
    lenient_mean_gapSize=mean_gapSize-1;
    trueSL=zeros(size(sl));
    for y=2:Rdim_seg_SL
        
        if(sl(y,1)==1) %White Col.
            height_count=height_count+1;
        end
        %         totalsSL=sum(sl,2,'omitnan');
        %         moving_total=0;
        %
        if(sl(y,1)==0)  %Black Col
            if(sl(y-1,1)==1)%Prev is white
                if(y>12) && (y<Rdim_seg_SL-12) && (y>round(height_count/2))%not near borders
                    for j=-9:9
                        candidate_rows=[candidate_rows sum(binarized_Segment((y-1-round(height_count/2))+j,:))];%binarized_Segment
                        candidate_locations=[candidate_locations (y-1-round(height_count/2))+j];
                    end
                    
                    for r=1:9
                        toRemove=[];
                        for s=1:size(candidate_rows,2)-1
                            if ((candidate_locations(1,s+1)-candidate_locations(1,s)<lenient_mean_gapSize))
                                if candidate_rows(1,s+1) > candidate_rows(1,s)
                                    toRemove = [toRemove; s];
                                else
                                    toRemove = [toRemove; s+1];
                                end
                            end
                        end
                        candidate_rows=candidate_rows(setdiff(1:length(candidate_rows),toRemove));
                        candidate_locations=candidate_locations(setdiff(1:length(candidate_locations),toRemove));
                    end
                    
                    [~,INDEX] = max(candidate_rows);
                    loc=candidate_locations(1,INDEX);
                    %sl(y-1:y-1-height_count,:)=0;
                    
                    toRemove=[];
                    for s=-7:7
                        if (loc>7)&&(loc<Rdim_seg_SL-7) && (trueSL(loc+s,1)==1)
                            toRemove = [toRemove; INDEX];
                            break;
                        end
                    end
                    candidate_rows=candidate_rows(setdiff(1:length(candidate_rows),toRemove));
                    candidate_locations=candidate_locations(setdiff(1:length(candidate_locations),toRemove));
                    
                    [M,INDEX] = max(candidate_rows);
                    loc=candidate_locations(1,INDEX);
                    trueSL(loc:loc,:)=1;
                    height_count=0;
                    candidate_rows=[];
                    candidate_locations=[];
                end
            end
        end
        
    end
    % ---------------------Getting average gap size-------------------------
    gap_count=1;
    height_count=0;
    totalSLCount=0;
    current_gap=[];
    SL_Locations = zeros(1,65);
    for y=2:Rdim_seg_SL
        if(trueSL(y,1)==1) %White Col.
            if(trueSL(y-1,1)==0) %Prev is white
                current_gap=[current_gap, gap_count];
                gap_count=0;
            end
            height_count=height_count+1;
        end
        
        if(trueSL(y,1)==0)  %Black Col
            if(trueSL(y-1,1)==1) %Prev is white
                totalSLCount=totalSLCount+1;
                SL_Locations(1,totalSLCount)=y-round(height_count/2);
                height_count=0;
            end
            gap_count=gap_count+1;
        end
    end
    mean_gapSize=median(current_gap);
    
    %#################################Data blocks##########################
    %.......... `.=. ,,,,,,,,,,,,,,
    % ---------------------Getting mean values of blocks between scale lines-------------------------
    gap_count=0;
    data_guess=(gray_Segment);
    for y=1:Rdim_seg
        if(trueSL(y,1)==0) && (y<Rdim_seg)
            gap_count=gap_count+1;
        end
        if(trueSL(y,1)==1) && (y<Rdim_seg)
            data_guess(y-1:y+1,:)=data_guess(y-1:y+1,:)-mean2(data_guess(y-1:y+1,:));
            gray_gap=data_guess(y-gap_count:y,1:Cdim_seg_SL);
            fun4=@(block_struct) mean2(block_struct.data)*ones(size(block_struct.data));
            data_gap =blockproc( gray_gap, [round(((gap_count)+1)/2),round(Cdim_seg_SL/40)], fun4);
            data_guess(y-gap_count:y,1:Cdim_seg_SL)=data_gap(:,:);
            gap_count=0;
        end
        if (y==Rdim_seg)
            data_guess(y-1:y,:)=data_guess(y-1:y,:)-mean2(data_guess(y-1:y,:));
            gray_gap=data_guess(y-gap_count:y,1:Cdim_seg_SL);
            fun4=@(block_struct) mean2(block_struct.data)*ones(size(block_struct.data));
            data_gap =blockproc( gray_gap, [round(((gap_count)+1)/2),round(Cdim_seg_SL/40)], fun4);
            data_guess(y-gap_count:y,1:Cdim_seg_SL)=data_gap(:,:);
            gap_count=0;
        end
    end
    data_thresh=data_guess;
    for y=1:Rdim_seg
        
        for x=1:Cdim_seg
            peak=double(round(0.75*(max(max(data_guess(:,x))))));
            if (data_thresh(y,x)<peak)
                data_thresh(y,x)=0;
            end
        end
    end
    
    regMax = imregionalmax(histeq(gray_SegmentUSER));
    cleared_regMax = bwareaopen(regMax, 50);
    
    mask = cleared_regMax;
    
    rp_data = regionprops(cleared_regMax, 'PixelIdxList', 'Orientation', 'Eccentricity');
    % Get high eccentricity and orientations at 90 and 0 degrees
    rp_data = rp_data(((abs([rp_data.Orientation]) > 88)|(abs([rp_data.Orientation]) < 2))& (abs([rp_data.Eccentricity]) > 0.9));
    
    mask(vertcat(rp_data.PixelIdxList)) = 0;
    data_only=mask;
    
    %     % ---------------------Thresholding-------------------------
    %     peak=double(round(0.50*(max(max(data_guess)))));
    %     data_thresh=im2bw(data_guess,double((peak)/255));
    %
    %     % ---------------------Finding largest CC-------------------------
    data_only_mask=imdilate(data_only,strel('octagon', 30));
    CC = bwconncomp(data_only_mask);
    numPixels = cellfun(@numel,CC.PixelIdxList);
    [biggest,idx] = max(numPixels);
    data_region=false(size(data_thresh));
    data_region(CC.PixelIdxList{idx}) = 1;
    
    data_only=data_region & data_only;
    %     %     % ---------------------Cropping Grayscale with CC-mask -------------------------
            data_only_grey=(gray_Segment);
            for y=1:Cdim_seg_SL
                for z=1:Rdim_seg_SL
                    if(data_only(z,y)==0)
                        data_only_grey(z,y)=0;
                    end
                end
            end

    % ---------------------Finding max pixel per col candidates, in masked grayscale-------------------------
    candidata_data_locations=[];
    candidata_data_totals=[];
    %data_only=imdilate(data_only,strel('octagon', 6));
    DATA_FINAL=zeros(size(data_only));
    for y=1:Cdim_seg
        for x=3:Rdim_seg-3
            moving_total=round((sum(data_only_grey(x-2,y)))+(sum(data_only_grey(x-1,y)))+(sum(data_only_grey(x,y)))+(sum(data_only_grey(x+1,y)))+(sum(data_only_grey(x+2,y)))/5);
            %if moving_total>0
            candidata_data_locations=[candidata_data_locations x];
            candidata_data_totals=[candidata_data_totals moving_total];
            %end
        end
        
        
        % ---------------------Placing max pixel per col, on blank image-------------------------
        [M,INDEX] = max(candidata_data_totals);
        if INDEX >0
            loc=candidata_data_locations(1,INDEX);
            if loc > 4
                DATA_FINAL(loc,y) = 1;
            end
        end
        candidata_data_totals=[];
        candidata_data_locations=[];
    end
    dataPoints = zeros(1,Cdim_seg_SL);
    for y=1:Cdim_seg_SL
        FOUND=0;
        for z=1:Rdim_seg_SL
            if(DATA_FINAL(z,y)==1) && (z>3)
                dataPoints(y)=Rdim_seg_SL-z;
                FOUND=1;
            end
        end
        if FOUND == 0
            dataPoints(y)=NaN;
        end
    end
    
    
    for y=6:Cdim_seg_SL-6
        moving_thresh=1.234*std(dataPoints(y-5:y+5),'omitnan');
        if(abs(dataPoints(y)-mean(dataPoints(y-5:y+5),'omitnan'))>moving_thresh)
            dataPoints(y) = NaN;
        end
    end

    
    %dataPoints_noNaN = inpaint_nans(dataPoints,3);
    %dataPoints_noNaN = (dataPoints);
    data_values = [data_values dataPoints];
    DIVISIONS=12;
    %figure, plot(dataPoints)
    %figure, plot(dataPoints_noNaN)
    DIVISIONS_deltaN_Values = [];
    avg_data=size(dataPoints);
    DIVISIONS_Interval=floor(size(dataPoints,2)/DIVISIONS);
    for j=0:DIVISIONS-1
        MEAN = median(dataPoints(1,(j*DIVISIONS_Interval+1):(j+1)*DIVISIONS_Interval),'omitnan');
        avg_data(1,(j*DIVISIONS_Interval+1):(j+1)*DIVISIONS_Interval)= MEAN;
    end
    avg_data_Values = [avg_data_Values avg_data];
    for j=2:DIVISIONS-1
        deltaN=(avg_data(1,floor((j-0.5)*(DIVISIONS_Interval)))-avg_data(1,floor((j-1.5)*(DIVISIONS_Interval)))) / DIVISIONS_Interval;
        DIVISIONS_deltaN_Values = [DIVISIONS_deltaN_Values deltaN];
    end
    
    deltaN_values = [deltaN_values DIVISIONS_deltaN_Values];
    if (i==3)
        %deltaN_back=(dataPoints_noNaN(1,2)-dataPoints_noNaN(size(dataPoints_noNaN,2))) / double(size(dataPoints_noNaN,2));
        background=mean(deltaN_values, 'omitnan');
    end
    if (i>11)
        for j=1:numel(DIVISIONS_deltaN_Values),
        percent_Increase = (1- DIVISIONS_deltaN_Values(j)/background)*100;
        percent_Increase_Values=[percent_Increase_Values percent_Increase];
        end
    end
    %plot(percent_Increase_Values)
    %
    % ---------------------Finding max pixel per col candidates, per white row-------------------------
    candidata_sl_locations=[];
    candidata_sl_totals=[];
    DATA_MASK=imdilate(DATA_FINAL,strel('octagon', 6));
    % ---------------------Cropping Grayscale with CC-mask -------------------------
    dedata_gray_Segment=gray_Segment;
    SL=zeros(size(gray_Segment));
    for y=1:Rdim_seg
        if trueSL(y,13)==1 && y>9 && y<Rdim_seg-9 %White Col.
            for x=1:Cdim_seg
                for z=-8:8
                    moving_total=round((sum(dedata_gray_Segment(y+z-1,x)))+(sum(dedata_gray_Segment(y+z,x)))+(sum(dedata_gray_Segment(y+z+1,x)))/3);
                    if moving_total>0
                        candidata_sl_locations=[candidata_sl_locations y+z];
                        candidata_sl_totals=[candidata_sl_totals moving_total];
                    end
                end
                % ---------------------Placing max pixel per col, on blank image-------------------------
                [M,INDEX] = max(candidata_sl_totals);
                if INDEX >0
                    loc=candidata_sl_locations(1,INDEX);
                    SL(loc,x) = 1;
                end
                candidata_sl_totals=[];
                candidata_sl_locations=[];
            end
        end
        
    end
    SL = bwareaopen(SL, 6);
    DATA_MASK=imdilate(DATA_FINAL,strel('octagon', 3));
    % ---------------------Cropping Grayscale with CC-mask -------------------------
    dedata_Segment=SL;
    SL=SL & ~DATA_MASK;
    trueSL=imdilate(SL,strel('octagon', 6));
    SL=SL & trueSL;

    % ---------------------Rebuild Segments-------------------------
    if(i==1)
        idealSL=SL;
        input_ideal=input_Segment;
        binarized_ideal=binarized_Segment;
        desprocketed_original=desprocketed_Segment;
        data=DATA_FINAL;
    else
        idealSL=cat(2,idealSL,SL);
        input_ideal=cat(2,input_ideal,input_Segment);
        binarized_ideal=cat(2,binarized_ideal,binarized_Segment);
        desprocketed_original=cat(2,desprocketed_original,desprocketed_Segment);
        data=cat(2,data,DATA_FINAL);
    end
    
    segment_message=['Segment_',num2str(i),'/',num2str(size(hmMiddle,2)-1),' processed.'];
    disp(segment_message)
end

%imwrite(idealHM,strcat(path,'04b-Hour Markers.bmp'),'bmp');

dimension_ideal=size(gray);
Rdim_ideal=dimension_ideal(1);
Cdim_ideal=dimension_ideal(2);

%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%       _        _        _        _        _        _               _        _        _        _        _        _
%      / /\     / /\     / /\     / /\     / /\     / /\           /\ \     /\ \     /\ \     /\ \     /\ \     /\ \
%     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \         /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \
%    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \       / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \
%   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/       \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \
%  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \     Ideal     / / /     \ \ \   / / /     \ \ \   / / /     \ \ \
%  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / /  _  Image  _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /
%   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\       /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /
%    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /       \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /
%     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /         \  / /   \  / /   \  / /   \  / /   \  / /   \  / /
%      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/           \/_/     \/_/     \/_/     \/_/     \/_/     \/_/
%
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
ideal=zeros(size(inputBW));
for x=1:Rdim_ideal
    for y=1:Cdim_ideal
        if (idealSprockets(x,y)==1)
            ideal(x,y)=1;
        end
        if (idealSL(x,y)==1)
            ideal(x,y)=1;
        end
        if (hm(x,y)==1)
            ideal(x,y)=1;
        end
        if (data(x,y)==1)
            ideal(x,y)=1;
        end
    end
end

imwrite(ideal,strcat(path,'05-Ideal.bmp'),'bmp');
imwrite(desprocketed_original,strcat(path,'desprocketed.bmp'),'bmp');

ideal_fat=imdilate(ideal,strel('octagon', 3));
ideal_comp=imcomplement(ideal_fat);
seed=zeros(size(input));
seed=double(input)/255;
for x=1:Rdim_ideal-1
    for y=1:Cdim_ideal-1
        if (ideal_comp(x,y)==1)
            seed(x,y,:)=NaN;
        end
    end
end

synthetic = zeros(size(input));
synthetic(:,:,1) = inpaint_nans(seed(:,:,1),2);
synthetic(:,:,2) = inpaint_nans(seed(:,:,2),2);
synthetic(:,:,3) = inpaint_nans(seed(:,:,3),2);

data_fat=imdilate(data,strel('octagon', 3));
blurDATA = 3 * imgaussfilt(data_fat,5);

for x=1:Rdim_ideal
    for y=1:Cdim_ideal
        if (data(x,y)~=0)
             synthetic(x,y,1)=synthetic(x,y,1)+blurDATA(x,y);
             synthetic(x,y,2)=synthetic(x,y,2)+blurDATA(x,y);
             synthetic(x,y,3)=synthetic(x,y,3)+blurDATA(x,y);
         end
    end
end

% imwrite(seed,strcat(path,'05-seed.bmp'),'bmp');
imwrite(synthetic,strcat(path,'05-synthetic.bmp'),'bmp');

%FACT: I need more coffee c[_]
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%       _        _        _        _        _        _               _        _        _        _        _        _
%      / /\     / /\     / /\     / /\     / /\     / /\           /\ \     /\ \     /\ \     /\ \     /\ \     /\ \
%     / /  \   / /  \   / /  \   / /  \   / /  \   / /  \         /  \ \   /  \ \   /  \ \   /  \ \   /  \ \   /  \ \
%    / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \ / / /\ \       / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \ / /\ \ \
%   / / /\_\// / /\_\// / /\_\// / /\_\// / /\_\// / /\_\/       \/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \\/_/\ \ \
%  / / /     \ \ \   / / /     \ \ \   / / /     \ \ \     Over-     / / /     \ \ \   / / /     \ \ \   / / /     \ \ \
%  \ \ \   _ / / /  _\ \ \   _ / / /  _\ \ \   _ / / /  _  lay    _  \ \ \ _   / / /_  \ \ \ _   / / /_  \ \ \ _   / / /
%   \ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\\ \ \/ /\       /\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / //\ \/ / /
%    \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ / \ \ \/ /       \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / / \ \/ / /
%     \ \  /   \ \  /   \ \  /   \ \  /   \ \  /   \ \  /         \  / /   \  / /   \  / /   \  / /   \  / /   \  / /
%      \_\/     \_\/     \_\/     \_\/     \_\/     \_\/           \/_/     \/_/     \/_/     \/_/     \/_/     \/_/
%
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�
overlay=desprocketed_original;
for x=1:Rdim_ideal
    for y=1:Cdim_ideal
        if (idealSL(x,y)==1)
            overlay(x,y,1)=255;
        end
        if (hm(x,y)==1)
            overlay(x,y,3)=255;
        end
        if (data(x,y)==1)
            overlay(x,y,2)=255;
        end
    end
end
imwrite(overlay,strcat(path,'07-Overlay.bmp'),'bmp');
imwrite(data,strcat(path,'07-Data.bmp'),'bmp');

%                                       __-------___
%                    �             _-~~             ~~-_  �����,��,�`����,��,�`����,�,��,��,�`����,��,�
%                    � *         _-~                    /~-_
%             /`\_/`\ �        /~  \                   /    \  �����,��,�`����,��,�`����,�,��,��,�`����,��,�
%           /|  ?|| ?|        /      \_______________/        \
%          | |___||__|      /       /     ///\\\     \          \  �����,��,�`����,��,�`����,�,��,��,�`����,��,�
%          |          \    /      /       |#| |#|      \          \
%          |    (==\/==) /______/          \\\ "        \_________ \  �����,��,�`����,��,�`����,�,��,��,�`����,��,�
%          |          / /        \        _  \\\       /            \
%           \         \^\         \      |#|  |#|     /               \  �����,��,�`����,��,�`����,�,��,��,�`����,��,�
%             \    &    ||          \_____\\\///____/      _-_       //\__/
%               \       ||------_-~~-_ ------------- \ --/~   ~\    || __/
%                 ~-----||====/~     |##################|       |/~~~~~
%                  (_(__/  ./     /                    \_\      \.
%                         (_(___/                         \_____)_)

disp('Operation complete.')
%                         ???????????  END  -??????????