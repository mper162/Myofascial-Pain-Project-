clc; clear; close all;

% --- Load mask ---
mask_file = 'traps_R.nii.gz';
mask = niftiread(mask_file);
info = niftiinfo(mask_file);

% --- Ensure binary ---
mask = mask > 0;

% --- Define 3D structuring element ---
se = strel('disk', 1);   % radius = 1 voxel, 3D erosion

% --- Apply erosion to entire 3D volume ---
mask_eroded = imerode(mask, se);

% --- Save eroded mask ---
output_file_uncompressed = 'traps_R_eroded.nii';  % save as .nii first
niftiwrite(uint8(mask_eroded), output_file_uncompressed);

% --- Compress to .nii.gz ---
% gzip(output_file_uncompressed);  % creates traps_eroded.nii.gz
% delete(output_file_uncompressed); % remove uncompressed file
