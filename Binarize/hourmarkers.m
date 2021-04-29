function output=hourmarkers(input,HOUR_MARKER_COLOUR)
%°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º°`°º¤ø,¸,ø¤°º¤ø,¸¸,ø¤º°`°º¤ø,¸¸,ø¤º° 
%Function to find hour markers on historic cosmic ray recordings, using
%   average intensity and Hough transform.
%Format:
%               output=hourmarkers(input)
%
%Input:
%               input: Unprocessed historic cosmic ray recording.
%
%Output:
%               output: An image containing only the ideal hour markers.
%
%Usage example:
%               idealHM = hourmarkers(binarized);

%Scale lines
%--------------------------------------------------------------------------
path='F:\Ground Truth\Generate\Binarize\';
dimension=size(input);
Rdim=dimension(1);
Cdim=dimension(2);
ideal_HM=zeros(size(input));

if(HOUR_MARKER_COLOUR==1)
    totalsHM=sum(input,1,'omitnan');
    for x=1:Cdim,
        if (totalsHM(x) > (mean(totalsHM)))%0.75
            ideal_HM(:, x) = 1;
        end
    end
end
if(HOUR_MARKER_COLOUR==0)
    totalsHM=sum(input,1,'omitnan');
    for x=1:Cdim-5,
        moving_total=round((totalsHM(x)+totalsHM(x+1)+totalsHM(x+2)+totalsHM(x+3)+totalsHM(x+4)+totalsHM(x+5))/5);
        average=mean(totalsHM);
        if (moving_total < (1.0*average))
            ideal_HM(:,x:x+4) = 1;
        end
    end
end
% 
% %   ??????????? Skeletonization--------------------------
% ideal_HM=imdilate(ideal_HM,strel('rectangle', [3 3]));
% ideal_HM=bwmorph(ideal_HM,'skel',Inf);
% ideal_HM=bwmorph(ideal_HM, 'spur', 10);
% 
% totals=sum(ideal_HM,1,'omitnan');
% thinned_ideal_HM=false(size(input));
% for x=1:numel(totals)-1,
%     if (totals(x) > (0.5 *max(totals)))
%         thinned_ideal_HM(:, x) = 1;
%     end
% end
% 

current_HM=[];
hmCandidates= zeros(1,14);
hmLocations= zeros(1,14);
totalHMCount=0;
width_count=0;
processed_ideal_HM=ideal_HM;
%processed_ideal_HM=imdilate(ideal_HM,strel('rectangle', [5 5]));
for x=2:Cdim,
    if(processed_ideal_HM(1,x)==1) %White Col.
        current_HM = [current_HM; sum(input(:,x),'omitnan');];
        width_count=width_count+1;
    end
    
    if(processed_ideal_HM(1,x)==0)  %Black Col
        if(processed_ideal_HM(1,x-1)==1) %Prev is white
            totalHMCount=totalHMCount+1;
            hmLocations(1,totalHMCount)=x-round(width_count/2);
            hmCandidates(1,totalHMCount)=sum(current_HM)/width_count;%-round(width_count/2);
            current_HM=[];
            width_count=0;
        end
    end
end

toRemove=[];
for i=1:size(hmLocations,2)
    if (hmLocations(1,i) ==0)
        toRemove = [toRemove; i];
    end
end
hmCandidateSUMS=hmCandidates(setdiff(1:length(hmCandidates),toRemove));
hmLocations=hmLocations(setdiff(1:length(hmLocations),toRemove));

[M,INDEX] = min(hmCandidateSUMS);
HM=zeros(size(input));
loc=hmLocations(1,INDEX);
if(loc+5<Cdim-5)
HM(:,loc:loc+5) = 1;
else
HM(:,loc:loc+1) = 1;
end
%   ??????????? Finding hour markers via Hough--------------------------
% start_angle = -2.5; % >-90 && <90
% end_angle = 2.5; % >-90 && <90
% theta_resolution = 0.75;
% [H T R] = hough(input,'RhoResolution',0.9, 'Theta', start_angle:theta_resolution:end_angle);
% P  = houghpeaks(H,1,'threshold',ceil(max(H(:))));
% lines = houghlines(input, T, R, P, 'MinLength',3);
%
% for k=1:length(lines);
%     %---------------Removing Hour Markers via Bresenham Algorithm-------------
%     [x, y] = bresenham(lines(k).point1(1), lines(k).point1(2), ...
%         lines(k).point2(1), lines(k).point2(2));
%     ideal_HM(y, x) = 1;
% end

output=HM;
end