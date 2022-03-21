% FUNCTION NAME:
%   generateH
%
% DESCRIPTION:
%   Generates an H-matrix for a computed tomography imaging spectrometer
%   (CTIS) system
%
% INPUT:
%   x           - Horizontal spatial dimension (# of columns) of input hyperspectral cube
%   y           - Vertical spatial dimension (# of rows) of input hyperspectral cube
%   z           - Spectral dimension of input hyperspectral cube (= # of spectral
%                   channels)
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
%
% OUTPUT:
%   H - Sparse system matrix describing the CTIS system
%
% ASSUMPTIONS AND LIMITATIONS:
%
%
% REVISION HISTORY
%   21/3/2022 - Mads Svanborg Peters
%       * Initial implementation
%

function H = generateH(x,y,z,b1,b2,shift,allOrders,diff_sens,illum,sigma_psf)

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

z = shift*z;
g_x = 3*x+2*z+2*(b1+b2-shift);      % # of rows in simulated CTIS image g
g_y = 3*y+2*z+2*(b1+b2-shift);      % # of columns in simulated CTIS image g

%rows = g_x*g_y;                    % # of rows in H
cols = x*y*z;                       % # of columns in H

Hc = cell(1,cols);                  % Allocate memory for H

iter = 1;
iter_k = 1;

if allOrders
    % For all 9 diffraction spots
    for k = 1:shift:z
        tic
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
        
        row_offset = zeros(1,9);
        col_offset = zeros(1,9);
        for j = 1:y
            for i = 1:x
                % Set row and column offsets for diffraction spots
                row_offset(:) = diff_offsets.pos(:,1) + i;
                col_offset(:) = diff_offsets.pos(:,2) + j;
                % Calculate factors for projected voxels
                factors = illum(iter_k).*diff_sens(:,iter_k);
                
                if bool_psf
                    % Generate sparse diffraction image with PSF
                    diff_image = zeros(g_x,g_y,'double');           % Allocate memory
                    for l = 1:9
                        diff_image(row_offset(l),col_offset(l)) = factors(l);
                    end
                    diff_image = imgaussfilt(diff_image,sigma_psf); % Apply PSF
                    diff_image = sparse(diff_image);
                else
                    % Generate sparse diffraction image without PSF
                    diff_image = sparse(row_offset,col_offset,factors',g_x,g_y);
                end
                % Insert vectorized sparse diff_image as column in H
                Hc{iter} = diff_image(:);
                
                iter = iter + 1;
            end
        end
        iter_k = iter_k + 1;
        disp(['Completion percentage = ',num2str(round(100.*k/z)),' %, ','Estimated time remaining: ',num2str(round(toc*(z-k))),' s'])
    end
    disp('Completed construction of H.')
    H = [Hc{:}]; % Output H
else
    % For 5 diffraction spots
    for k = 1:shift:z
        tic
        % Indices of the 5 diffraction orders
        diff_offsets.pos(1,1:2) = [(g_x-x)/2 - x-b1-(k-1), (g_y-y)/2];  % Top
        
        diff_offsets.pos(2,1:2) = [(g_x-x)/2, (g_y-y)/2 - y-b1-(k-1) ]; % Left
        diff_offsets.pos(3,1:2) = [(g_x-x)/2, (g_y-y)/2];               % Middle
        diff_offsets.pos(4,1:2) = [(g_x-x)/2, (g_y-y)/2 + y+b1+(k-1)];  % Right
        
        diff_offsets.pos(5,1:2) = [(g_x-x)/2 + x+b1+(k-1), (g_y-y)/2];  % bottom
        
        row_offset = zeros(1,5);
        col_offset = zeros(1,5);
        for j = 1:y
            for i = 1:x
                % Set row and column offsets for diffraction spots
                row_offset(:) = diff_offsets.pos(:,1) + i;
                col_offset(:) = diff_offsets.pos(:,2) + j;
                % Calculate factors for projected voxels
                factors = illum(iter_k).*diff_sens(:,iter_k);
                
                if bool_psf
                    % Generate sparse diffraction image with PSF
                    diff_image = zeros(g_x,g_y,'double');           % Allocate memory
                    for l = 1:9
                        diff_image(row_offset(l),col_offset(l)) = factors(l);
                    end
                    diff_image = imgaussfilt(diff_image,sigma_psf); % Apply PSF
                    diff_image = sparse(diff_image);
                else
                    % Generate sparse diffraction image without PSF
                    diff_image = sparse(row_offset,col_offset,factors',g_x,g_y);
                end
                % Insert vectorized sparse diff_image as column in H
                Hc{iter} = sparse(diff_image(:));
                iter = iter + 1;
            end
        end
        iter_k = iter_k + 1;
        disp(['Completion percentage = ',num2str(round(100.*k/z)),' %, ','Estimated time remaining: ',num2str(round(toc*(z-k))),' s'])
    end
    disp('Completed construction of H.')
    
    H = [Hc{:}]; % Output H
end

end

