
scanProtocol = struct('Ns', 300, 'NTheta', 360, 'mmPerSample', 1);

%% 1a
ellipse = struct('h', 0, 'k', 0,'a', 10, 'b', 10, 'alpha', 0);
sinogramA = parallelBeamProjection(ellipse, scanProtocol);
figure; imshow(sinogramA', []); colorbar;
export_fig('hw2p1a','-pdf');

%% 1b 
ellipse = struct('h', 80, 'k', 0,'a', 10, 'b', 10, 'alpha', 0);
sinogramB = parallelBeamProjection(ellipse, scanProtocol);
figure; imshow(sinogramB', []); colorbar;
export_fig('hw2p1b','-pdf');
%% 1c
ellipse = struct('h', 80, 'k', 0,'a', 10, 'b', 10, 'alpha', 0);
sinogramC = parallelBeamProjection(ellipse, scanProtocol);
figure; imshow(sinogramC', []); colorbar;
export_fig('hw2p1c','-pdf');
%% 1d
ellipse = struct('h', 5, 'k', 0,'a', 10, 'b', 10, 'alpha', 0);
sinogramD = parallelBeamProjection(ellipse, scanProtocol);
figure; imshow(sinogramD', []); colorbar;
export_fig('hw2p1d','-pdf');
%% 1e
muSkull = 1;  % [cm-1]
muBrain = 0.2; % [cm-1]
ellipseSkull = struct('h', 0, 'k', 0,'a', 70, 'b', 100, 'alpha', 0);
ellipseNegSkull = struct('h', 0, 'k', 0,'a', 65, 'b', 95, 'alpha', 0);
ellipseBrain = ellipseNegSkull;

%% Projections
sinogramSkull = parallelBeamProjection(ellipseSkull, scanProtocol);
sinogramNegSkull = parallelBeamProjection(ellipseNegSkull, scanProtocol);
sinogramBrain = parallelBeamProjection(ellipseBrain, scanProtocol);
%%
sinogramBrainPhantomE = muSkull *(sinogramSkull - sinogramNegSkull)+ ...
     muBrain * sinogramBrain;
% im = iradon(sinogramBrainPhantom', linspace(1e-5*180/pi,(2*pi-1e-5)*180/pi, NTheta));

figure; imshow(sinogramBrainPhantomE', []); colorbar;
export_fig('hw2p1e','-pdf');
