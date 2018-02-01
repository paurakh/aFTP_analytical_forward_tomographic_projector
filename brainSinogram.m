%% ellipse
a = 80; % [mm]
b = 100; % [mm]
h = 0; % [mm]
k = 0; % [mm]
alpha = 0; %[radians]

rotatingDetSource;
sinogramBrain = 1 * sinogram;

%%
a = 75; % [mm]
b = 95; % [mm]
h = 0; % [mm]
k = 0; % [mm]
alpha = 0; %[radians]

rotatingDetSource;
sinogramBrainNeg = 1 * sinogram;

%%

a = 75; % [mm]
b = 95; % [mm]
h = 0; % [mm]
k = 0; % [mm]
alpha = 0; %[radians]

rotatingDetSource;
sinogramWater = 0.2* sinogram;

sinogramFinal = sinogramBrain - sinogramBrainNeg + sinogramWater;
%%
figure; imshow(sinogramFinal', [0 2]); colorbar
im = iradon(sinogramFinal', theta*180/pi);
figure; imshow(im, []); colorbar
