clc; clear; close all;

%% ================== User Inputs ==================
fat_signal_path   = 'fat_ideal_neck.nii.gz';
water_signal_path = 'water_ideal_neck.nii.gz';

mask_files = {
    'traps_eroded_l.nii';
    'traps_eroded_R.nii'
};

out_excel = 'dixon_traps_signal_quantification.xlsx';

%% ================== Load Maps ==================
fat_signal   = double(niftiread(fat_signal_path));     
water_signal = double(niftiread(water_signal_path));   

%% ================== Output Table ==================
results = table( ...
    'Size',[numel(mask_files) 3], ...
    'VariableTypes',{'string','double','double'}, ...
    'VariableNames',{'MaskName','Fat_pct','Water_pct'} ...
);

%% ================== Loop Over Masks ==================
for i = 1:numel(mask_files)

    maskName = mask_files{i};
    mask = niftiread(maskName) > 0;

    % ---------- Extract signals ----------
    fat_vals   = fat_signal(mask);
    water_vals = water_signal(mask);

    % Remove invalid values
    valid_idx = isfinite(fat_vals) & isfinite(water_vals) & (fat_vals >= 0) & (water_vals >= 0);

    fat_vals   = fat_vals(valid_idx);
    water_vals = water_vals(valid_idx);

    % ---------- Calculate percentages ----------
    total_signal = fat_vals + water_vals;
    fat_pct   = mean(fat_vals ./ total_signal * 100);
    water_pct = mean(water_vals ./ total_signal * 100);

    % ---------- Store ----------
    results.MaskName(i)  = string(erase(maskName,'.nii.gz'));
    results.Fat_pct(i)   = fat_pct;
    results.Water_pct(i) = water_pct;

end

%% ================== Write Excel ==================
writetable(results, out_excel);