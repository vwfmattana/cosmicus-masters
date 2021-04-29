function output=sprockets(input)
%����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,�`����,�,��,��,�`����,��,� 
%Function to find sprockets on historic cosmic ray recordings, using
%correlation.
%Format:
%               output=sprockets(input)
%
%Input:
%               input: Unprocessed historic cosmic ray recording.
%
%Output:
%               output: An image containing only the ideal sprockets.
%
%Usage example:
%               idealSprockets = sprockets(inputHSV);

%Sprockets - Finding via Correlation
%--------------------------------------------------------------------------
path='/home/vincent/cosmicus/code/Binarize';
dimension=size(input);
inputHSV=rgb2hsv(double(input));
Rdim=dimension(1);%Mdim
Cdim=dimension(2);%Ndim
CORRELATION_THRESHOLD=0.9;
idealSprocketTemplate=imread(strcat(strcat(path,'/'),'idealTemplate.bmp'));
thresholdedIdealSprocketTemp = im2bw(idealSprocketTemplate(:,:,1), 0.5);

sprocketTemplate=imread(strcat(strcat(path,'/'),'template0.bmp'));
thresholdedSprocketTemp = im2bw(sprocketTemplate(:,:,1), 0.5);
correlationOutput = normxcorr2(sprocketTemplate(:,:,1), inputHSV(:,:,3));
thresholdedCorrelation = im2bw(correlationOutput(:,:),(CORRELATION_THRESHOLD*max(correlationOutput(:)))); % adjust

%   ??????????? Closing & Clustering---------------------------
sprocketPoints = imclose(thresholdedCorrelation, strel('line',4, 90));
[labels, num] = bwlabel(sprocketPoints, 8);

idealSprocketLocations = false(size(sprocketPoints));

for i = 1:num
   [r, c] = find(labels==i); 
   x = round(mean(c));
   y = round(mean(r));
   idealSprocketLocations(y,x) = 1;
end

%   ??????????? Finding Ideal Sprocket---------------------------
for m=1:Cdim,
    for n=1:Rdim,
        if (idealSprocketLocations(n,m)>0 && n>(size(thresholdedSprocketTemp,1)) && m>(size(thresholdedSprocketTemp,2)));
            startRow = n-(size(thresholdedSprocketTemp,1));
            startCol = m-(size(thresholdedSprocketTemp,2));
            realTemplate = input(startRow:startRow+(size(thresholdedSprocketTemp,1)-1),startCol:startCol+(size(thresholdedSprocketTemp,2)-1),:);
        end
    end
end
imwrite(realTemplate,strcat(path,'TRUETEMP.bmp'),'bmp');
realSprocketTemplateHSV=rgb2hsv(double(realTemplate));
thresholdedRealTemplate = im2bw(realTemplate, 0.75);

CORRELATION_THRESHOLD=0.70;
correlationOutputFinal = normxcorr2(realSprocketTemplateHSV(:,:,3), inputHSV(:,:,3));
thresholdedCorrelationFinal = im2bw(correlationOutputFinal(:,:),(CORRELATION_THRESHOLD*max(correlationOutputFinal(:)))); % adjust

%   ??????????? Erosion----------------------------------------
sprocketPointsFinal=imerode(thresholdedCorrelationFinal,strel('diamond',1));

%   ??????????? Closing & Clustering---------------------------
sprocketPoints = imclose(sprocketPointsFinal, strel('line',4, 90));
[labels, num] = bwlabel(sprocketPoints, 8);

idealSprocketLocationsFinal = false(size(sprocketPoints));

for i = 1:num
   [r, c] = find(labels==i); 
   x = round(mean(c));
   y = round(mean(r));
   idealSprocketLocationsFinal(y,x) = 1;
end

%   ??????????? Placing Ideal Sprockets---------------------------
idealSprockets = false(size(idealSprocketLocationsFinal));
for m=1:Cdim,
    for n=1:Rdim,
        if (idealSprocketLocationsFinal(n,m)>0 && n>(size(thresholdedRealTemplate,1)) && m>(size(thresholdedRealTemplate,2))+4 && ((n<80)||(n>Rdim-80)));
            startRow = n-(size(thresholdedRealTemplate,1));
            startCol = m-(size(thresholdedRealTemplate,2))-5;
            idealSprockets(startRow:startRow+(size(thresholdedIdealSprocketTemp,1)-1),startCol:startCol+(size(thresholdedIdealSprocketTemp,2)-1))= thresholdedIdealSprocketTemp;
        end
    end
end

output=idealSprockets;
end