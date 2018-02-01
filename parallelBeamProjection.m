function [ sinogram ] = parallelBeamProjector( ellipse, scanProtocol )
% Parallel beam forward Projector
% Paurakh Rajbhandary
%%
Ns = scanProtocol.Ns;        % number of detectors
NTheta = scanProtocol.NTheta;    % number of views
mmPerSample = scanProtocol.mmPerSample; % 1 mmPerSample of detector

% define reference x-coordinate array
x = linspace(-Ns/2* mmPerSample, Ns/2*mmPerSample, Ns);


yDet = -1500 * ones(1,Ns);                 % reference detector y-coord array
ySource = 1500 * ones(1,Ns);               % reference source y-coord array
theta = linspace(1e-5,2*pi-1e-5, NTheta);  % rotation angle array

% catenate reference source and detector array
xySource = [x' ySource'];
xyDet = [x' yDet'];

% figure;
% Loop for every view
for iTheta = 1:length(theta)
    
    % rotation matrix
    rotMatrix = [cos(theta(iTheta)) -sin(theta(iTheta));
        sin(theta(iTheta)) cos(theta(iTheta))];
    
    % rotate both source and detector array
    xyDetRot = xySource * rotMatrix;
    xySourceRot = xyDet * rotMatrix;
    
    % calculate x and y intercept for these rotated array
    q = (xyDetRot(:,2) ./ xyDetRot(:,1) - xySourceRot(:,2) ./ xySourceRot(:,1)) ...
        ./ (1./xyDetRot(:,1) - 1 ./xySourceRot(:,1));
    p = 1 ./ (1 ./xyDetRot(:,1)- xyDetRot(:,2) ./ (xyDetRot(:,1) .*q));
    
    % compute coordinates of intersection of x-rays and object
    [x1 x2 y1 y2 ] = lineEllipseMatrix(ellipse.a, ellipse.b, ellipse.h, ellipse.k, ellipse.alpha, p, q);
    
    % set non-intersecting lines (with comple value) to be zero
    % so that the distance between them is also zero
    cond1 = (imag(x1) ~= 0) | (imag(x2) ~= 0);
    x1(cond1) = 0;
    y1(cond1) = 0;
    x2(cond1) = 0;
    y2(cond1) = 0;
    % compute the distance
    sinogram(iTheta,:) = sqrt((x1-x2).^2+(y1-y2).^2)';
end


end

