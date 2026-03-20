
clc; clear all; 


input_folder = 'UTE_data';          % folder containing UTE NIfTIs
output_file  = 'merged_ute.nii.gz'; % output merged file
tes_file     = 'TEs_2.txt';           % output TE list

% Call function
merge_UTE_by_TE(input_folder, output_file, tes_file);


function merge_UTE_by_TE(input_folder, output_file, tes_file)

    % --- Get all .nii and .nii.gz files ---
    files = dir(fullfile(input_folder, '*.nii*'));
    if isempty(files)
        error('No NIfTI files found in folder: %s', input_folder);
    end

    file_list = fullfile(input_folder, {files.name});
    nFiles = numel(file_list);
    TEs = zeros(nFiles,1);

    % --- Extract TE from metadata ---
    for i = 1:nFiles
        fname = file_list{i};
        if endsWith(fname, '.gz')
            gunzip(fname);
            fname = erase(fname, '.gz');
        end

        info = niftiinfo(fname);

        % Try EchoTime first
        if isfield(info, 'EchoTime') && ~isempty(info.EchoTime)
            TEs(i) = info.EchoTime;  
        else
            % fallback: parse from Description
            desc = info.Description;
            te_val = regexp(desc, 'TE\s*[:=]\s*([0-9.]+)', 'tokens');
            if ~isempty(te_val)
                TEs(i) = str2double(te_val{1}{1});
            else
                error('Could not extract TE from metadata for file: %s', fname);
            end
        end

        % Update with unzipped filename
        file_list{i} = fname;
    end

    % --- Sort by TE ---
    [TEs_sorted, sort_idx] = sort(TEs);
    file_list = file_list(sort_idx);

    % --- Load first volume ---
    ref_vol  = niftiread(file_list{1});
    ref_info = niftiinfo(file_list{1});

    % Initialize 4D array
    merged_data = zeros([size(ref_vol), nFiles], class(ref_vol));

    % --- Load all volumes ---
    for i = 1:nFiles
        merged_data(:,:,:,i) = niftiread(file_list{i});
    end

    % --- Update header for 4D ---
    ref_info.ImageSize = size(merged_data);
    ref_info.PixelDimensions = [ref_info.PixelDimensions 1];
    ref_info.Datatype = class(merged_data);

    % --- Write 4D NIfTI ---
    niftiwrite(merged_data, output_file, ref_info, 'Compressed', true);

    % --- Write TE list ---
    fid = fopen(tes_file, 'w');
    fprintf(fid, '%g ', TEs_sorted);
    fclose(fid);

    fprintf('✅ Merged %d UTE datasets into %s (sorted by TE)\n', nFiles, output_file);
    fprintf('📝 TE list written to %s\n', tes_file);
end
