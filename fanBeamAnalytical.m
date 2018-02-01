function [ sinogram ] = fanBeamAnalytical(shape, SDD, SID, fOV, nDet, nViews, detectorOffset, rotationRange, nSSource, nSDet)

%% define fanbeam
% SDD = 100; % [cm]  SDD > SID + FOV/2
% SID = 60; % [cm]  == D (yes its true)
% fOV = 50 ; % [cm]

fanAngle = atan(fOV/2/SID)*180/pi * 2; % [degrees]
% rotationRange =[0 360]; % [90-fanAngle/2 360]; % [degrees]


nSView = 1;

% detectorOffset = 0;  % [fraction] ex: 1/4
focalSpotSize = 1/10 ; %[cm]
% nViews = 2000;
% nDet = 888;


nDet0 = nDet;
nViews0 = nViews;
nViews = nViews * nSView;
nDet = nDet * nSDet;


focalSpotAngle = focalSpotSize / SID * 180/pi;

rotationAngle = linspace(rotationRange(1), rotationRange(2), nViews+1); rotationAngle = rotationAngle(1:end-1);

rotationAngle1 = zeros(nSSource, nViews);
for i = 1:size(rotationAngle,2)
     a = linspace(rotationAngle(i)-focalSpotAngle/2, rotationAngle(i)+focalSpotAngle/2, nSSource);
     rotationAngle1(:,i) = a;
end
rotationAngle2 = reshape(rotationAngle1, [1 size(rotationAngle1,1) * size(rotationAngle1,2)]);
rotationAngle2(find(rotationAngle2<0)) = rotationAngle2(find(rotationAngle2<0)) + 360;

rotationAngle(find(rotationAngle<0)) = rotationAngle(find(rotationAngle<0)) + 360;
%thetaRanage = [alpha - (fanAngle/2) alpha - (fanAngle/2)];
xSource = -SID*cos(rotationAngle2*pi/180);
ySource = -SID*sin(rotationAngle2*pi/180);
relFanAngle = linspace(-fanAngle/2, +fanAngle/2, nDet);
% quarter Detector Offset
offsetAngle = detectorOffset*(relFanAngle(2)-relFanAngle(1));
relFanAngle = relFanAngle + offsetAngle;

xDetRangeRel = (SDD)*cos(relFanAngle*pi/180)-SID;
yDetRangeRel = (SDD)*sin(relFanAngle*pi/180);

XsYsXdYd = zeros(4,nViews*nDet*nSSource);

xSourceGrid = reshape(repmat(xSource,[nDet 1]), [1 nDet*nViews*nSSource]);
ySourceGrid = reshape(repmat(ySource,[nDet 1]), [1 nDet*nViews*nSSource]);
XsYsXdYd(1,:) = xSourceGrid;

XsYsXdYd(2,:) = ySourceGrid;


rotationAngleGrid = repmat(rotationAngle',[1 nSSource nDet]);


xDetRangeRelGrid = repmat(permute(xDetRangeRel,[1 3 2]),[nViews nSSource 1]);
yDetRangeRelGrid = repmat(permute(yDetRangeRel,[1 3 2]),[nViews nSSource 1]);

xDetRangeAbsGrid = cos((rotationAngleGrid)*pi/180) .* xDetRangeRelGrid ...
    - sin((rotationAngleGrid) *pi/180) .* yDetRangeRelGrid;
yDetRangeAbsGrid = sin((rotationAngleGrid)*pi/180) .* xDetRangeRelGrid ...
    + cos((rotationAngleGrid) *pi/180) .* yDetRangeRelGrid;

xDetRangeAbsGrid = reshape(permute(xDetRangeAbsGrid, [ 3 2 1]), [1 nDet*nViews*nSSource]);
yDetRangeAbsGrid = reshape(permute(yDetRangeAbsGrid, [ 3 2 1]), [1 nDet*nViews*nSSource]);

xyDet = [xDetRangeAbsGrid; yDetRangeAbsGrid];

XsYsXdYd(3,:) = xDetRangeAbsGrid; %reshape(xDetRangeAbsGrid', [1 size(xDetRangeAbsGrid,1)*size(xDetRangeAbsGrid,2)]);
XsYsXdYd(4,:) = yDetRangeAbsGrid; %reshape(yDetRangeAbsGrid', [1 size(yDetRangeAbsGrid,1)*size(yDetRangeAbsGrid,2)]);

clear xDetRangeAbsGrid;
clear yDetRangeAbsGrid
clear xSourceGrid;
clear ySourceGrid;

p = ((XsYsXdYd(2,:)-XsYsXdYd(4,:)) ./ (XsYsXdYd(1,:)-XsYsXdYd(3,:)) .* XsYsXdYd(1,:) - XsYsXdYd(2,:)) ...
    ./ ((XsYsXdYd(2,:)-XsYsXdYd(4,:)) ./ (XsYsXdYd(1,:)-XsYsXdYd(3,:)));

q = -p .* (XsYsXdYd(2,:) - XsYsXdYd(4,:)) ./ (XsYsXdYd(1,:)-XsYsXdYd(3,:));

x = []; y = [];
for i= 1:length(shape)
    [x1 x2 y1 y2] = lineEllipseMatrix(shape(i).a, shape(i).b, shape(i).h, shape(i).k, shape(i).alpha, p, q );
    if (not(isnan(x1)))
        x = [x; x1 ; x2];
        y = [y; y1 ; y2];
    end
end
xSGrid = repmat(squeeze(XsYsXdYd(1,:)), [size(x,1) 1]);
ySGrid = repmat(squeeze(XsYsXdYd(2,:)), [size(y,1) 1]);

x(imag(x) ~= 0) = xSGrid(imag(x) ~= 0);
y(imag(y) ~= 0) = ySGrid(imag(y) ~= 0);


%%
d = zeros(size(x,1)/2,size(x,2));

d(1,:) = sqrt((x(1,:)-x(2,:)) .^2 + (y(1,:) - y(2,:)) .^2);
for iD = 2: (size(x,1)/2)
    
    d(iD,:) = sqrt((x(iD*2-1,:)-x(iD*2,:)) .^2 + (y(iD*2-1,:) - y(iD*2,:)) .^2);
    d(1,:) = d(1,:) - d(iD,:);
end

for iShape = 1:size(shape,2)
    d(iShape,:) = d(iShape,:) * shape(iShape).mu;
end

d = sum(d,1);
% 
% d = zeros(1, nDet*nViews*nSSource);
% %%
% % if (not(isempty(x)))
%     x1 = x;
%     y1 = y;
%     [x y] = sortFromSource(x, y ,squeeze(XsYsXdYd(1,:)), squeeze(XsYsXdYd(2,:)));  
%     
%     %%
%     for i = 1:(size(x,1)-1)
%         xMid = mean([x(i,:) ; x(i+1,:)]);
%         yMid = mean([y(i,:) ; y(i+1,:)]);
%         mu = getMaterial4SegmentMatrix(xMid, yMid, shape);
%         d = d + mu .* sqrt((x(i, :)-x(i+1, :)) .^ 2 + (y(i, :)-y(i+1, :)) .^ 2);
%     end
% % end
% %sinogram = reshape(d, [nDet nViews]);

sinogram = reshape(d, [nDet nSSource nViews]);
sinogram = squeeze(mean(sinogram,2));
sinogram = reshape(sinogram, [nSDet nDet0 nViews]);
sinogram = squeeze(mean(sinogram,1));

% nn = 3;
% nViews = nn*nViews;
% sinogram1 = zeros(nDet, nViews); 
% for i = 1:nn
% sinogram1(:,i:nn:end) = sinogram;
% end
% sinogram = sinogram1;
% sinogram = reshape(sinogram, [nDet nViews0 nSView]);
% sinogram = squeeze(mean(sinogram,3));

nViews = nViews0;
nDet = nDet0;
%figure; imagesc( rotationAngle,relFanAngle,sinogram)

end

