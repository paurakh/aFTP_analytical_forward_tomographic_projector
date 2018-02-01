
%% ellipse
a = 3;
b = 1;
h = 2;
k = 2;
alpha = pi/4;
%%
Ns = 300;
NTheta = 360;
cmPerSample = 0.1; % 1 mmPerSample

x = linspace(-Ns/2* cmPerSample, Ns/2*cmPerSample, Ns);

yDet = -150 * ones(1,Ns);
ySource = 150 * ones(1,Ns);
theta = linspace(1e-5,2*pi-1e-5, NTheta);
xySource = [x' ySource'];
xyDet = [x' yDet'];
figure;
p = 0; q = inf;
for iTheta = 1:length(theta)
    
    rotMatrix = [cos(theta(iTheta)) -sin(theta(iTheta));
                sin(theta(iTheta)) cos(theta(iTheta))];

    xyDetRot = xySource * rotMatrix;
    xySourceRot = xyDet * rotMatrix;
    
    q = (xyDetRot(:,2) ./ xyDetRot(:,1) - xySourceRot(:,2) ./ xySourceRot(:,1)) ...
        ./ (1./xyDetRot(:,1) - 1 ./xySourceRot(:,1));
    p = 1 ./ (1 ./xyDetRot(:,1)- xyDetRot(:,2) ./ (xyDetRot(:,1) .*q));
    [x1 x2 y1 y2 ] = lineEllipseMatrix(a, b, h, k, alpha, p, q)
    
    cond1 = (imag(x1) ~= 0) | (imag(x2) ~= 0);
    x1(cond1) = 0; y1(cond1) = 0; x2(cond1) = 0; y2(cond1) = 0;
    
    sinogram(iTheta,:) = sqrt((x1-x2).^2+(y1-y2).^2)';
    
end

figure; imshow(sinogram', [0 2]); colorbar
im = iradon(sinogram', theta*180/pi);
figure; imshow(im, [])

% 
% x1, y1, x2, y2
% 
% x1 / p + y1 / q = 1
% 1 p' + y1/x1 q' = 1/x1
% 
% 1 p' + y2/x2 q' = 1/x2
% 
% q' (y1/x1 - y2/x2) = 1/x1 - 1/x2
% q' = (1/x1 - 1/x2) /(y1/x1 - y2/x2)
% 
% q = (y1/x1 - y2/x2)/ (1/x1 - 1/x2)
% p = 1/ (1/x1- y1/(x1*q));