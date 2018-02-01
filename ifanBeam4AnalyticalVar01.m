function [ im ] = ifanBeam4AnalyticalVar01(sinogram, SDD, SID, fOV, fOVR, xSize, ySize, nDet, nViews, detectorOffset, rotationRange )
%% back projection

% step 1
% RBetaI' (n*alpha) = RBetaI (n*alpha) * D * cos( n * alpha)
fanAngle = atan(fOV/2/SID)*180/pi * 2; % [degrees]


relFanAngle = linspace(-fanAngle/2, +fanAngle/2, nDet);
offsetAngle = detectorOffset*(relFanAngle(2)-relFanAngle(1));
relFanAngle = relFanAngle + offsetAngle;

rotationAngle = linspace(rotationRange(1), rotationRange(2), nViews+1); rotationAngle = rotationAngle(1:end-1);
rotationAngle(find(rotationAngle<0)) = rotationAngle(find(rotationAngle<0)) + 360;

xSource = -SID*cos(rotationAngle*pi/180);
ySource = -SID*sin(rotationAngle*pi/180);

sinogram1 = sinogram* (SID ^ 2) .* (repmat(cos(relFanAngle*pi/180)',[1 nViews]) .^2);
% sinogram1 = sinogram; %* (SID^2) .* (repmat(cos(relFanAngle*pi/180)',[1 nViews]) .^2);

% step 2 (using FFT)
% QBetaI (n*alpha) = RBetaI' (n*alpha) ** g (n * alpha)
n = (-2*nDet):(2*nDet);
alpha = abs(relFanAngle(2)-relFanAngle(1))*pi/180;
a = (1/(8*alpha^2)) * ( n ==0 );
b = (0) * (mod(n,2) == 0);
c = -(1/2)*((1./(pi*sin(n*alpha))).^2) .* (mod(n,2) == 1);
c (find(isnan(c))) = 0 ;
g = a .* not(isnan(a)) + b .* not(isnan(b))+ c .* not(isnan(c));

G = fftshift(fft(g,512))';

ha = (1/(4*alpha^2)) * ( n ==0 );
hb = (0) * (mod(n,2) == 0);
hc = (-1./(n.^2 * pi^2 * alpha^2)) .* (mod(n,2) == 1);
hc (find(isnan(hc))) = 0 ;
h = ha + hb + hc;
H = fftshift(fft(h, 512))';

% xSize = 512/512;
% ySize = 512/512;

im = zeros(1, xSize*ySize);
% deltaBeta = abs(rotationAngle(2)-rotationAngle(1))*pi/180;
xCenter =  0;
yCenter = 0;
xCorL = linspace((xCenter-fOVR/2),(xCenter+fOVR/2), xSize);
yCorL = linspace((yCenter-fOVR/2),(yCenter+fOVR/2), ySize);
[xCorGrid yCorGrid] = meshgrid(xCorL,yCorL);
clear xCorL
clear yCorL;

xCorGrid = xCorGrid(:)';
yCorGrid = yCorGrid(:)';
r1 =  zeros(2 * nDet-1 - size(sinogram,1),1);
r1 = [zeros(2^nextpow2(size(sinogram,1)) - size(sinogram,1),1); r1];

hanning = hann(length(g));
FKernel = fftshift((fft(g))) .* (hanning');

kernel = ifft(ifftshift(FKernel));
% kernel(:) = 1;
%kernel = g;
% kernel = kernel .^ 2;
testq = zeros(size(sinogram,1),nViews);

for iView = 1:nViews
    xSourceSel = xSource(iView);
    ySourceSel = ySource(iView);
    %     r = [r1 ; sinogram1(:,iView)];
    
    %     q = conv(r, kernel);
    %     q = q(length(r) + length(sinogram(:,1)) +1 : length(r) + 2* length(sinogram(:,1)));
    %
    q = conv(squeeze(sinogram1(:,iView)), kernel .^2, 'same')';
%      q = sinogram1(:,iView)';
    
    %     testq(:,iView) = q';
    
    % q1 = conv(r, g);
    % q1 = q1(length(r) + length(sinogram(:,1)) +1 : length(r) + 2* length(sinogram(:,1)));
    
    absFanAngle = linspace(rotationAngle(iView)-fanAngle/2, rotationAngle(iView)+fanAngle/2, nDet);
    absFanAngle(find(absFanAngle<0)) = absFanAngle(find(absFanAngle<0)) + 360;
    
    atanAngle = atan(abs((yCorGrid - ySourceSel) ./ (xCorGrid-xSourceSel)))*180/pi;
    absAngle = atanAngle;
    ind = find( ((yCorGrid - ySourceSel) > 0) & ( (xCorGrid - xSourceSel) <0 ) ) ;
    absAngle(ind) = 180 - absAngle(ind);
    ind =  find( ((yCorGrid - ySourceSel) < 0) & ( (xCorGrid - xSourceSel) <0 ));
    absAngle(ind) = 180 + absAngle(ind);
    ind = find( ((yCorGrid - ySourceSel) < 0) & ((xCorGrid - xSourceSel) - 0) > 0);
    absAngle(ind) = 360 - absAngle(ind);
    clear atanAngle;
%     absAngle = gather(absAngle);
    if (absFanAngle(1) >= (360-fanAngle))
        absFanAngle(absFanAngle <= fanAngle) = absFanAngle(absFanAngle <= fanAngle) + 360;
        absAngle(absAngle<= fanAngle) = absAngle(absAngle <= fanAngle) + 360;
    end
    
    index = (absAngle-absFanAngle(1))/(absFanAngle(end)-absFanAngle(1));
    index((index<=0) | (index>=1)) = 0;
    index = index * (nDet-1) +1;
    i1 = floor(index);
    i2 = ceil(index);
    
    val1 = q(i1);
    val2 = q(i2);
    angle1 = absFanAngle(i1);
    angle2 = absFanAngle(i2);
    if ((angle2 - angle1) ~= 0)
        term1 = (absAngle -angle1)./(angle2 - angle1);
    else
        term1 = 1;
    end
%     term1 = (absAngle -angle1)./(angle2 - angle1);
    val = (1-term1) .^ 2 .* val1 + term1 .^2 .* val2;
    %          val = val1;
    %     val = (1 - ((absAngle -angle1)./(angle2 - angle1)) ) .^2 .* val1; ...
    %         + ((absAngle -angle1)./(angle2 - angle1)) .^2 .* val2;
    %     val = val1 + ((absAngle -angle1)./(angle2 - angle1)).^2 .*((val2 + val1));
    
    %     val = val1;
    %      val = val2;
    c1 = (absAngle >= angle1) & (absAngle <= angle2);
    c2 = (absAngle >= angle2) & (absAngle <= angle1);
    c = (c1 | c2);
    val = val .* c;
    val(isnan(val)) = 0;
    L = sqrt((yCorGrid-ySourceSel).^2+(xCorGrid-xSourceSel).^2);
    %     im = im + val;
    
    im = im + (1./L.^2) .^ 2 .* val;
    if (mod(iView,10) == 0)
        dispProgress(iView,nViews);
        %figure(4); imagesc(flipud(reshape(im,[xSize ySize])));
    end
    
end
% backprojection
im = reshape(im', [sqrt(size(im,2)) sqrt(size(im,2))]);
im = (flipud(im* (1.6*pi/(nDet*nViews)) .^ 2));
%figure(4); imshow(im, [-0.1/100, 0.1/100]);
%figure(6); imshow(im, [0.2-0.1/100,0.2+ 0.1/100]);

% disp(im(end/2,end/2));
% val = im(end/2,end/2);
% k = val;
% disp(k);
% disp(0.2/k);
% k1 = k/(nViews*nDet);
% disp(k1);
end

