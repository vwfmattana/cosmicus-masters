function output=scalelines(inputBW)
%°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°
%Function to find scale lines on historic cosmic ray recordings, using
%average intensity and Hough transform.
%Format:
%               output=scalelines(inputBW)
%
%Input:
%               input: Unprocessed historic cosmic ray recording.
%
%Output:
%               output: An image containing only the ideal scale lines.
%
%Usage example:
%               idealSL = scalelines(binarized);

%Scale lines
%--------------------------------------------------------------------------
path='F:\Ground Truth\Generate\Binarize\';
dimension=size(inputBW);
Rdim=dimension(1);
Cdim=dimension(2);
ideal_SL=zeros(size(inputBW));
%figure, imshow(inputBW)
%   ??????????? Getting sum of Scale Lines---------------------------
totalsSL=sum(inputBW,2,'omitnan');
for x=1:Rdim-2,
    moving_total=(totalsSL(x)+totalsSL(x+1)+totalsSL(x+2))/3;
    average=mean(totalsSL);
    roof=max(totalsSL);
    if (moving_total > (average)) %|| (moving_total > (average))
        ideal_SL(x:x+2,:) = 1;
    end
end
% %   ??????????? Finding Scale Lines via Hough--------------------------
% start_angle = -90; % >-90 && <90
% end_angle = -89; % >-90 && <90
% theta_resolution = 0.7;
% 
% [H T R] = hough(inputBW,'RhoResolution',0.9, 'Theta', start_angle:theta_resolution:end_angle);
% P  = houghpeaks(H,500,'threshold',ceil(0.001*max(H(:))),'NHoodSize',[NHoodS 2]);
% lines = houghlines(inputBW, T, R, P, 'MinLength',3);
% 
% %   ??????????? Writing buffed Lines via Bresenham Algorithm-------------
% for l=1:length(lines);
%     [x, y] = bresenham(lines(l).point1(1), lines(l).point1(2), ...
%         lines(l).point2(1), lines(l).point2(2));
%     ideal_SL(y, x) = 1;
% end

% ideal_SL=imdilate(ideal_SL,strel('rectangle', [3 3]));
% ideal_SL=bwmorph(ideal_SL,'skel',Inf);
% ideal_SL=bwmorph(ideal_SL, 'spur', 10);

current_SL=[];
slCandidates= zeros(65,1);
slLocations= zeros(65,1);
totalSLCount=0;
height_count=0;
%processed_ideal_SL=imdilate(ideal_SL,strel('rectangle', [2 2]));
processed_ideal_SL=ideal_SL;
for x=2:Rdim,
    if(processed_ideal_SL(x,1)==1) %White Row.
        current_SL = [current_SL; sum(inputBW(x,:),'omitnan');];
        height_count=height_count+1;
    end
    
    if(processed_ideal_SL(x,1)==0)  %Black Row
        if(processed_ideal_SL(x-1,1)==1) %Prev is white
            totalSLCount=totalSLCount+1;
            slLocations(totalSLCount,1)=x-round(height_count/2);
            slCandidates(totalSLCount,1)=sum(current_SL)/height_count;%-round(width_count/2);
            current_SL=[];
            height_count=0;
        end
    end
end

toRemove=[];
for i=1:size(slLocations,1)
    if (slLocations(i,1) ==0)
        toRemove = [toRemove; i];
    end
end
slCandidates=slCandidates(setdiff(1:length(slCandidates),toRemove));
slLocations=slLocations(setdiff(1:length(slLocations),toRemove));

[M,INDEX] = max(slCandidates);
SL=false(size(inputBW));
loc=slLocations(INDEX,1);
SL(loc,:) = 1;

output=SL;
end

