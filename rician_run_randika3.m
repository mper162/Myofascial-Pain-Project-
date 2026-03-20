
clc; clear all;
delete(gcp('nocreate'));   % prevent parallel pool errors

%% --- Load Inputs ---
% Load magnitude data

%nii_mag = load_untouch_nii('ute_mag.nii.gz');

 % min TE = 0.028 ms
% imData.mag = double(nii_mag.img);  % Ensure it's double
% imData.mag = squeeze(imData.mag); % (x, y, z, echoes)
imData.images = squeeze(double(niftiread('merged_ute.nii.gz')));

% Load echo times
TE = load('TEs_1.txt');  % Assume in ms
TE = TE ./ 1000; %%  ms to s very important!!!
imData.TE = TE(:)';    % Row vector
imData.TE = imData.TE;    %Randy addition 

% Load ROI mask
%nii_mask = load_untouch_nii('VM_t1_dseg.nii.gz');
nii_mask_1 = niftiread('traps_L_eroded.nii');
nii_mask = nii_mask_1(:,:,:,1);
%nii_mask = squeeze(double(niftiread('VM_t1_dseg.nii.gz')));  % randy
roi.mask = logical(nii_mask);
 % figure;
 % imshow(roi.mask(:,:,21))
% Metadata
imData.FieldStrength = 3.0; % Tesla
%imData.SignalModel = 'rician'; % Optional field if needed by class

%%%%%%%%%%%%%%%*******************************************************

%%%%testing code

% Assuming 'imData.images' and 'roi.mask' are already loaded
% and contain your 4D image data and 3D logical mask, respectively.

% slice_number = 21; % Choose the slice you want to display
% te_index = 1;      % Choose the echo time for the main image (e.g., first TE)
% 
% % 1. Extract the main MRI image slice
% main_image_slice = imData.images(:, :, slice_number, te_index);
% 
% % 2. Extract the corresponding mask slice
% mask_slice = roi.mask(:, :, slice_number);
% 
% % 3. Create a new figure for the display
% figure;
% 
% % 4. Display the main MRI image slice
% %    Use '[]' to auto-scale the intensity values for optimal contrast.
% imshow(main_image_slice, []);
% colormap gray; % Set the colormap to grayscale for anatomical image
% 
% % 5. Hold the current plot to overlay the mask
% hold on;
% 
% % 6. Create an RGB overlay for the mask
% %    Where the mask is true (1), it will be green; otherwise, transparent.
% %    You can change the color by modifying the [R G B] values (e.g., [1 0 0] for red).
% colored_mask = zeros([size(mask_slice), 3]);
% colored_mask(mask_slice, 2) = 1; % Set the green channel to 1 where mask is true
% 
% % 7. Display the colored mask with transparency
% %    'AlphaData' controls transparency: 0 is fully transparent, 1 is fully opaque.
% %    Adjust the multiplier (e.g., 0.3) to control the mask's visibility.
% h_overlay = imshow(colored_mask);
% set(h_overlay, 'AlphaData', mask_slice * 0.3); % Mask becomes 30% opaque where true
% 
% % 8. Add a title for clarity
% title(['Slice ' num2str(slice_number) ' with VM Mask Overlay (TE ' num2str(imData.TE(te_index)*1000) 'ms)']);
% 
% % 9. Release the plot hold
% hold off;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%********************************************



%% --- Initialize Fitting Class ---

fp = FittingProcessorRician_3Comp();%randy

%% --- Run Fitting ---

maps = fp.MultistepFitImage(imData, roi);%% --- Save Results ---%randy



% %Save fat fraction as NIfTI
% nii_FF = nii_mag;
% nii_FF.img = single(maps.FF);
% save_untouch_nii(nii_FF, 'FF_map.nii.gz');
% 
% % Save individual component maps if needed
% nii_mag.img = single(maps.Water);
% save_untouch_nii(nii_mag, 'Water_map.nii.gz');
% 
% nii_mag.img = single(maps.Fat);
% save_untouch_nii(nii_mag, 'Fat_map.nii.gz');
% 
% disp('Fitting complete. Results saved.');

%%%%%%%%%%%%%%  new Save Format %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% --- Load Inputs (existing part of your script) ---
% imData.images = squeeze(double(niftiread('ute_mag.nii.gz')));
% TE = load('TEs.txt');  % Assume in ms
% imData.TE = TE(:)';    % Row vector
% nii_mask = niftiread('VM_t1_dseg.nii.gz');
% roi.mask = logical(nii_mask);
% imData.FieldStrength = 3.0; % Tesla
% % ... (rest of your script, including calling MultistepFitImage) ...

%% --- Save Results using niftiwrite ---

% Assuming 'nii_info' is already defined as in the previous example

% Save Fat Fraction map
niftiwrite(single(maps.FF), 'FF_map.nii.gz');

% Save Total Water map (as before)
niftiwrite(single(maps.Water), 'Water_map.nii.gz');

% Save Total Fat map (as before)
niftiwrite(single(maps.Fat), 'Fat_map.nii.gz');

% --- Add these lines to save individual water components ---

% Save Long-T2* Water Initial Signal (S0) map
niftiwrite(single(maps.Water_long), 'S0water_long_map.nii.gz');

% Save Short-T2* Water Initial Signal (S0) map
niftiwrite(single(maps.Water_short), 'S0water_short_map.nii.gz');

% You might also want to save the R2* maps for these components
%niftiwrite(single(maps.R2water_long), 'R2water_long_map.nii.gz');
%niftiwrite(single(maps.R2water_short), 'R2water_short_map.nii.gz');
niftiwrite(single(maps.Wsfraction), 'Wsfraction_map.nii.gz');
niftiwrite(single(maps.T2str_water_short), 'T2str_water_short_map.nii.gz');
niftiwrite(single(maps.T2str_water_long), 'T2str_water_long_map.nii.gz');


disp('Fitting complete. All specified results saved.');
%disp('Fitting complete. Results saved.');