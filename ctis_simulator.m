% FUNCTION NAME:
%   ctis_simulator
%
% DESCRIPTION:
%   Generates a simulated CTIS image based on a hyperspectral cube
%
% INPUT:
%   cube        - 3D hyperspectral cube with dimensions (x, y, z)
%   b1          - [Optional] Border between zeroth- and first-orders in # of pixels
%   b2          - [Optional] Border between first-orders and outer edge of CTIS image in # of pixels
%   shift       - [Optional] Pixel shift between projected spectral bands in the first-order
%                   diffraction spot
%   allOrders   - [Optional] Boolean parameter to show 9 diffraction order spots (true)
%                   or 5 diffraction order spots (false)
%   diff_sens   - [Optional] Wavelength depedent diffraction sensitivity,
%                   i.e. diffraction efficiency including transmission and sensor response
%                   of system
%   illum       - [Optional] Spectrum of the illuminant used during the
%                   acquistion of CTIS images
%   sigma_psf   - [Optional] Standard deviation of the applied point spread function
%                   (PSF)
%   noise       - [Optional] Standard deviation of the applied zero mean
%                   white Gaussian noise
%
% OUTPUT:
%   g - Simulated CTIS image
%
% ASSUMPTIONS AND LIMITATIONS:
%
%
% REVISION HISTORY
%   21/3/2022 - Mads Svanborg Peters
%       * Initial implementation
%
function g = ctis_simulator(cube,b1,b2,shift,allOrders,diff_sens,illum,sigma_psf,noise)

[x, y, z] = size(cube);             % Determine dimensions of cube
% Check inputs and assign default values where necessary 

if ~exist('b1','var') || isempty(b1)
    b1 = 1;                         % Default value
end

if ~exist('b2','var') || isempty(b2)
    b2 = 0;                         % Default value
end

if ~exist('shift','var') || isempty(shift)
    shift = 1;                      % Default value
end

if ~exist('allOrders','var') || isempty(allOrders)
    allOrders = false;              % Default value
end

if ~exist('diff_sens','var') || isempty(diff_sens)
    if allOrders
        diff_sens = ones(9,z);      % Default value for allOrders = 1
    else
        diff_sens = ones(5,z);      % Default value for allOrders = 0
    end
end

if ~exist('illum','var') || isempty(illum)
    illum = ones(1,z);              % Default value
end

bool_psf = true;                    % Boolean for PSF
if ~exist('sigma_psf','var') || isempty(sigma_psf)
    bool_psf = false;               % Default value
end

bool_noise = true;                  % Boolean for noise
if ~exist('noise','var') || isempty(noise)
    bool_noise = false;             % Default value
end

z = shift*z;
g_x = 3*x+2*z+2*(b1+b2-shift);      % # of rows in simulated CTIS image g
g_y = 3*y+2*z+2*(b1+b2-shift);      % # of columns in simulated CTIS image g

iter_k = 1;

diff_image = zeros(g_x,g_y,class(cube)); % Allocate memory for simulated CTIS image

if allOrders
    % 9 diffraction order spots
    for k = 1:shift:z
        
        % Indices of the 9 diffraction orders
        diff_offsets.pos(1,1:2) = [(g_x-x)/2 - x-b1-(k-1), (g_y-y)/2 - y-b1-(k-1) ];    % Top-left
        diff_offsets.pos(2,1:2) = [(g_x-x)/2 - x-b1-(k-1), (g_y-y)/2];                  % Top
        diff_offsets.pos(3,1:2) = [(g_x-x)/2 - x-b1-(k-1), (g_y-y)/2 + y+b1+(k-1) ];    % Top-right
        
        diff_offsets.pos(4,1:2) = [(g_x-x)/2, (g_y-y)/2 - y-b1-(k-1) ];                 % Left
        diff_offsets.pos(5,1:2) = [(g_x-x)/2, (g_y-y)/2];                               % Middle
        diff_offsets.pos(6,1:2) = [(g_x-x)/2, (g_y-y)/2 + y+b1+(k-1)];                  % Right
        
        diff_offsets.pos(7,1:2) = [(g_x-x)/2 + x+b1+(k-1), (g_y-y)/2 - y-b1-(k-1) ];    % Bottom-left
        diff_offsets.pos(8,1:2) = [(g_x-x)/2 + x+b1+(k-1), (g_y-y)/2];                  % bottom
        diff_offsets.pos(9,1:2) = [(g_x-x)/2 + x+b1+(k-1), (g_y-y)/2 + y+b1+(k-1) ];    % Bottom-right
        
        row_off = zeros(1,9); % Allocate memory for offsets
        col_off = zeros(1,9);
        
        row_off(:) = diff_offsets.pos(:,1)+1;
        col_off(:) = diff_offsets.pos(:,2)+1;
        
        for ii = 1:9
            % Calculate factors for the 1 to 9 diffraction spots
            factor = illum(iter_k).*diff_sens(ii,iter_k);
            % Iteratively fill up the CTIS image with the 9 diffraction
            % spots
            I = factor.* cube(:,:,iter_k) + diff_image(row_off(ii):(row_off(ii)+x-1), col_off(ii):(col_off(ii)+y-1));
            % Insert into CTIS image, diff_image
            diff_image(row_off(ii):(row_off(ii)+x-1), col_off(ii):(col_off(ii)+y-1)) = I;
        end
        iter_k = iter_k + 1;
    end
    
    g = diff_image; % Simulated CTIS image
    
    % Apply PSF and noise
    if bool_psf
        g = imgaussfilt(diff_image,sigma_psf);
    end
    if bool_noise
        g = g + noise .* randn(g_x,g_y);    % Addition of zero mean white Gaussian noise
        g(g<0) = 0;                         % Set negative values to zero
    end
else
    % 5 diffraction order spots
    for k = 1:shift:z
        % Indices of the 5 diffraction orders
        diff_offsets.pos(1,1:2) = [(g_x-x)/2 - x-b1-(k-1), (g_y-y)/2];                  % Top
        
        diff_offsets.pos(2,1:2) = [(g_x-x)/2, (g_y-y)/2 - y-b1-(k-1) ];                 % Left
        diff_offsets.pos(3,1:2) = [(g_x-x)/2, (g_y-y)/2];                               % Middle
        diff_offsets.pos(4,1:2) = [(g_x-x)/2, (g_y-y)/2 + y+b1+(k-1)];                  % Right
        
        diff_offsets.pos(5,1:2) = [(g_x-x)/2 + x+b1+(k-1), (g_y-y)/2];                  % bottom
        
        row_off = zeros(1,5); % Allocate memory for offsets
        col_off = zeros(1,5);
        
        row_off(:) = diff_offsets.pos(:,1)+1;
        col_off(:) = diff_offsets.pos(:,2)+1;
        
        for ii = 1:5
            % Calculate factors for the 1 to 9 diffraction spots
            factor = illum(iter_k).*diff_sens(ii,iter_k);
            % Iteratively fill up the CTIS image with the 9 diffraction
            % spots
            I = factor.* cube(:,:,iter_k) + diff_image(row_off(ii):(row_off(ii)+x-1), col_off(ii):(col_off(ii)+y-1));
            % Insert into CTIS image, diff_image
            diff_image(row_off(ii):(row_off(ii)+x-1), col_off(ii):(col_off(ii)+y-1)) = I;
        end
        iter_k = iter_k + 1;
    end
    
    g = diff_image; % Simulated CTIS image
    
    % Apply PSF and noise
    if bool_psf
        g = gaussfilt(diff_image,sigma_psf);
    end
    if bool_noise
        g = g + noise .* randn(g_x,g_y);    % Addition of zero mean white Gaussian noise
        g(g<0) = 0;                         % Set negative values to zero
    end
    
end

end
