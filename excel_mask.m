clc; clear; close all;

%% ================== User Inputs ==================
fat_frac_path    = 'FF_map.nii.gz';                  % 0–100 scale
short_frac_path  = 'Wsfraction_map.nii.gz';    % 0–100 scale

t2_short_path = 'T2str_water_short_map.nii.gz';
t2_long_path  = 'T2str_water_long_map.nii.gz';

mask_files = {
    'traps_L_eroded.nii'
     'traps_R_eroded.nii'
};

out_excel = 'traps_quantification_eroded.xlsx';

%% ================== Load Maps ==================
fat_frac     = double(niftiread(fat_frac_path));     
short_frac   = double(niftiread(short_frac_path));   

t2_short_map = double(niftiread(t2_short_path));
t2_long_map  = double(niftiread(t2_long_path));

% Derived
water_frac = 100 - fat_frac;
long_frac  = 100 - short_frac;

%% ================== Output Table ==================
results = table( ...
    'Size',[numel(mask_files) 10], ...
    'VariableTypes',{ ...
        'string','double','double','double','double', ...
        'double','double','double','double','double'}, ...
    'VariableNames',{ ...
        'MaskName','Fat_pct','Water_pct', ...
        'Water_Long_pct','Water_Short_pct', ...
        'Final_FA_pct','Final_Water_Long_pct','Final_Water_Short_pct', ...
        'Water_Long_T2_ms','Water_Short_T2_ms'} ...
);

%% ================== Loop ==================
for i = 1:numel(mask_files)

    maskName = mask_files{i};
    mask = niftiread(maskName) > 0;

    % ---------- Fat / Water ----------
    fat_vals   = fat_frac(mask);
    water_vals = water_frac(mask);

    fat_vals   = fat_vals(isfinite(fat_vals) & fat_vals >= 0 & fat_vals <= 100);
    water_vals = water_vals(isfinite(water_vals) & water_vals >= 0 & water_vals <= 100);

    fat_pct   = mean(fat_vals);
    water_pct = mean(water_vals);

    % ---------- Short / Long ----------
    short_vals = short_frac(mask);
    long_vals  = long_frac(mask);

    short_vals = short_vals(isfinite(short_vals) & short_vals >= 0 & short_vals <= 100);
    long_vals  = long_vals(isfinite(long_vals)  & long_vals >= 0 & long_vals <= 100);

    short_pct = mean(short_vals);
    long_pct  = mean(long_vals);

    % ---------- T2 ----------
    t2s = t2_short_map(mask);
    t2l = t2_long_map(mask);

    valid_short = isfinite(t2s) & t2s > 0 & t2s < 100 & t2s ~= 65535;
    valid_long  = isfinite(t2l) & t2l > 0 & t2l < 200 & t2l ~= 65535;

    mean_t2_short = mean(t2s(valid_short));
    mean_t2_long  = mean(t2l(valid_long));

    % ---------- Combined ----------
    final_fa_pct          = fat_pct;
    final_water_long_pct  = (water_pct * long_pct) / 100;
    final_water_short_pct = (water_pct * short_pct) / 100;

    % ---------- Store ----------
    results.MaskName(i)              = string(erase(maskName,'.nii'));
    results.Fat_pct(i)               = fat_pct;
    results.Water_pct(i)             = water_pct;
    results.Water_Long_pct(i)        = long_pct;
    results.Water_Short_pct(i)       = short_pct;
    results.Final_FA_pct(i)          = final_fa_pct;
    results.Final_Water_Long_pct(i)  = final_water_long_pct;
    results.Final_Water_Short_pct(i) = final_water_short_pct;
    results.Water_Long_T2_ms(i)      = mean_t2_long;
    results.Water_Short_T2_ms(i)     = mean_t2_short;

end

%% ================== Write ==================
writetable(results, out_excel);