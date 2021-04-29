% =========================================================================
% Vincent Mattana 21128707 2016 MSc IT/Computer Science
% Hour Markers
% =========================================================================
clear;
[filename, pathname] = ...
     uigetfile({'*.bmp;*.jpg;',...
 'Picture Files (*.bmp,*.jpg)';},'Select source image for processing');
input=imread(strcat(strcat(pathname,'\'),filename) );

%--------------------Paths-------------------------------------------------
path='F:\LosTheRed\Computer Science\Masters\cosmicus\2016\Artefact\Ground Truth\Generate\Hour Markers\'; 


%// Create horizontal line structuring element
se = strel('line', 50, 90);

%// Dilate the image with this structuring element
out = imerode(input, se);
imshow(out)
 
%BW = imrotate(input,90);
%BW=input;
% [H,theta,rho] = hough(BW);
% 
% [R, xp] = radon(BW, 85:95);
% imagesc(0:180, xp, R)
% plot(R(:, 19))
% P  = houghpeaks(H,2);
% imshow(H,[],'XData',theta,'YData',rho,'InitialMagnification','fit');
% xlabel('\theta'), ylabel('\rho');
% axis on, axis normal, hold on;
% plot(theta(P(:,2)),rho(P(:,1)),'s','color','white');
% 
% 
% lines = houghlines(BW,theta,rho,P,'FillGap',5,'MinLength',500);
% 
% figure, imshow(BW), hold on
% max_len = 0;
% for k = 1:length(lines)
%    xy = [lines(k).point1; lines(k).point2];
%    plot(xy(:,1),xy(:,2),'LineWidth',2,'Color','green');
% 
%    % Plot beginnings and ends of lines
%    plot(xy(1,1),xy(1,2),'x','LineWidth',2,'Color','yellow');
%    plot(xy(2,1),xy(2,2),'x','LineWidth',2,'Color','red');
% 
%    % Determine the endpoints of the longest line segment
%    len = norm(lines(k).point1 - lines(k).point2);
%    if (len > max_len)
%       max_len = len;
%       xy_long = xy;
%    end
% end
% % highlight the longest line segment
% plot(xy_long(:,1),xy_long(:,2),'LineWidth',2,'Color','red');
