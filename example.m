% Example file for demonstration of system matrix H generation, simulation
% of CTIS image g and reconstruction of hyperspectral cube f using
% expectation maximization (EM)


% Define inputs (see generateH.m and projectg.m for definitions)
x = 100;
y = 100;
z = 25;
b1 = 5;
b2 = 0;
shift = 1; 
allOrders = true;
sigma_psf = 1;
noise = 0.1;
% Load illumination, diffraction sensitivity and wavelength axis 
load('sensitivity25.mat');
load('halogen25.mat');
load('wavelength25.mat');

diff_sens = sensitivity;
illum = halogen_z;
wave = wavelength;

%cube = ones(x,y,z); % Simple white cube

% Load example hyperspectral cube of a colorchecker
load('cube_HSI_colorchecker.mat');
cube = double(cube);

figure(1)
imagesc(sum(cube,3));  axis('equal'); axis([-inf inf -inf inf]); title('Summation of spectral bands');

% Crop cube to desired dimensions
cube = cube(50:50+y-1, 100:100+x-1,round(linspace(10,200,z)));


%% Visualize diffraction sensitivity and illumination
figure(2);

names = {'(1,1)','(2,1)','(3,1)','(1,2)','(2,2)','(3,2)','(1,3)','(2,3)','(3,3)'};

nexttile(1)
plot(wave,diff_sens','Linewidth',2); grid on;
xlabel('Wavelength [nm]'); ylabel('Relative sensitivity');
title('Diffraction sensitivity','interpreter','latex')
legend(names,'Location','Northwest','NumColumns',3)
axis([-inf inf -inf inf])

nexttile(2)
plot(wave,illum,'LineWidth',2); grid on;
xlabel('Wavelength [nm]'); ylabel('Normalized intensity'); title('Halogen lamp');
axis([-inf inf -inf inf]); 

%% Example for ideal CTIS system (no diff_sens, illum, PSF or noise)
% Generate system matrix H
H = generateH(x,y,z,b1,b2,shift,allOrders);

% Simulate CTIS image g using system matrix H  (g = Hf)
g_ideal = H*cube(:);
% Simulate CTIS image g using CTIS simulator  (g = Hf + n)
g = ctis_simulator(cube,b1,b2,shift,allOrders);

% Determine dimensions of g:
g_y = size(g,1);
g_x = size(g,2);
% Reshape ideal g into image
g_ideal = reshape(g_ideal,g_x,g_y);

figure (3)
nexttile(1)
imagesc(g_ideal); colorbar; axis('equal'); axis([-inf inf -inf inf]); title('g = Hf');

nexttile(2);
imagesc(g); colorbar; axis('equal'); axis([-inf inf -inf inf]); title('g = Hf + n');

% Reconstruct hyperspectral cube using EM
iter = 10;
f = em(H,g,iter);
% Reshape vectorized cube to 3D cube
cube_em = reshape(f,x,y,z);

% Visaulize reconstructed cube
figure (4)
montage(cube_em,'DisplayRange',[-inf inf]); 
title('Visualization of reconstructed spectral bands'); colorbar;


%% Example including diffraction sensitivity, illumination, PSF and noise
% Generate system matrix H (PSF results in significant decrease in
% sparsity and increase in computation time and memomory requirements)
H = generateH(x,y,z,b1,b2,shift,allOrders,diff_sens,illum,sigma_psf);

% Simulate CTIS image g using system matrix H  (g = Hf)
g_noNoise = H*cube(:);
% Simulate CTIS image g using CTIS simulator  (g = Hf + n)
g = ctis_simulator(cube,b1,b2,shift,allOrders,diff_sens,illum,sigma_psf,noise);

% Determine dimensions of g:
g_y = size(g,1);
g_x = size(g,2);
% Reshape g into image
g_noNoise = reshape(g_noNoise,g_x,g_y);

% Visualization of simulated CTIS images
figure (5)
nexttile(1)
imagesc(g_noNoise); colorbar; axis('equal'); axis([-inf inf -inf inf]); title('g = Hf');

nexttile(2);
imagesc(g); colorbar; axis('equal'); axis([-inf inf -inf inf]); title('g = Hf + n');

% Reconstruct hyperspectral cube using EM
iter = 10;
f = em(H,g,iter);
% Reshape vectorized cube to 3D cube
cube_em = reshape(f,x,y,z);

% Visaulize reconstructed cube
figure (6)
montage(cube_em,'DisplayRange',[-inf inf]); 
title('Visualization of reconstructed spectral bands'); colorbar; 

