classdef FittingProcessorRician_3Comp
    % FittingProcessorRician_3Comp
    %   This class implements multi‐step image fitting routines for water–fat
    %   separation using a three‐component model. In this model, the signal is
    %   assumed to be a sum of fat, water_long and water_short components – each with
    %   its own R2 decay (R2fat, R2water_long, and R2water_short). The initialization
    %   for the two‐point search forces an equal split between water_long and water_short,
    %   but once the fitting is done the fitted water components are used as obtained.
    %   For sigma estimation the usual two‐component model is used.
    
    methods (Static)
        function [maps, sigmaFromRoi, sigmaFromFit] = FitImageSigma(imData, roi)
            % FitImageSigma estimates sigma from image data via voxel-wise Rician fitting.
            % This method uses the standard 2-component model (as in our earlier class)
            % for sigma estimation.
            %
            % Syntax:
            %   [maps, sigmaFromRoi, sigmaFromFit] = FittingProcessorRician_3Comp.FitImageSigma(imData, roi)
            %
            % Inputs:
            %   imData - structure with fields:
            %            .images         4-D array (nx,ny,nz,nTE)
            %            .TE             1-by-n vector (echo times in seconds)
            %            .FieldStrength  Scalar (Tesla)
            %   roi    - structure with field:
            %            .mask Logical (nx-by-ny-by-nz) defining region of interest
            %
            % Outputs:
            %   maps         - structure with parameter maps (FF, R2water, R2fat, sigma, S0, SNR)
            %   sigmaFromRoi - sigma estimated from ROI statistics
            %   sigmaFromFit - sigma estimated from voxel-wise fitting (corrected)
            
            TE = 1000 * imData.TE;
            %TE =  imData.TE;
            TE = reshape(TE, [], 1);
            FieldStrength = imData.FieldStrength;
            
            mask = roi.mask;
            data = imData.images;
            [nx, ny, nz, nTE] = size(data);
            
            % Preallocate maps as column vectors
            FFrician = zeros(nx*ny*nz, 1);
            R2water_rician = zeros(nx*ny*nz, 1);
            R2fat_rician   = zeros(nx*ny*nz, 1);
            sigmaEstimates = zeros(nx*ny*nz, 1);
            s0Estimates    = zeros(nx*ny*nz, 1);
            
            mask_r = mask(:);
            data_r = reshape(data, nx*ny*nz, nTE);
            RoiIdxs = mask_r > 0 & all(data_r, 2);
            data_roi = data_r(RoiIdxs, :);
            nb_voxels = size(data_roi, 1);
            
            % Estimate sigma from ROI using standard deviation across echoes.
            stdMaskedIm = std(data_roi, [], 1);
            [~, TEindex] = min(abs(imData.TE - 0.0023));
            sigmaFromRoi = stdMaskedIm(TEindex);
            
            % For sigma estimation we use option 2 (2-component model)
            algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(data_roi, TE, sigmaFromRoi, 2);
            
            FFrician_roi = zeros(nb_voxels, 1);
            R2water_roi   = zeros(nb_voxels, 1);
            R2fat_roi     = zeros(nb_voxels, 1);
            sigmaEstimates_roi = zeros(nb_voxels, 1);
            s0Estimates_roi    = zeros(nb_voxels, 1);
            
            disp('Estimating sigma from ROI (2-component model)...');
            tstart = tic;
            parfor i = 1:nb_voxels
                ydata = data_roi(i,:);
                algoparams = algoparamsArray(i);
                outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting_WithSigma(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                % The 2-component model returns parameters: [F; W; R2; sigma]
                FFrician_roi(i) = outparams.F / (outparams.F + outparams.W);
                R2water_roi(i) = outparams.R2water;
                R2fat_roi(i)   = outparams.R2fat;
                sigmaEstimates_roi(i) = outparams.sig;
                s0Estimates_roi(i) = outparams.F + outparams.W;
            end
            tstop = toc(tstart);
            disp("Sigma estimated in: " + tstop + "s");
            
            FFrician(RoiIdxs,:) = FFrician_roi;
            R2water_rician(RoiIdxs,:) = R2water_roi;
            R2fat_rician(RoiIdxs,:) = R2fat_roi;
            sigmaEstimates(RoiIdxs,:) = sigmaEstimates_roi;
            s0Estimates(RoiIdxs,:) = s0Estimates_roi;
            
            maps.FFrician = reshape(FFrician, nx, ny, nz);
            maps.R2water   = reshape(R2water_rician, nx, ny, nz);
            maps.R2fat     = reshape(R2fat_rician, nx, ny, nz);
            maps.sigma     = reshape(sigmaEstimates, nx, ny, nz);
            maps.S0        = reshape(s0Estimates, nx, ny, nz);
            maps.SNR       = s0Estimates ./ sigmaEstimates;
            
            validSigma = maps.sigma(maps.sigma > 0);
            sigmaFromFit = mean(validSigma);
            correctionData = 1.163;
            sigmaFromFit = sigmaFromFit * correctionData;
        end
        
        function maps = FitImage(imData, sigma)
            % FitImage performs voxel-wise fitting on a given slice using a fixed sigma.
            % This routine estimates the three-component model:
            %   [F; Wl; Ws; R2fat; R2water_long; R2water_short]
            % where F is the fat signal and the water signal is split into water_long (Wl)
            % and water_short (Ws). The initialization for the 2-point search forces an equal
            % (half and half) split between Wl and Ws. Once fitting is done, the fitted values
            % are returned as computed.
            %
            % Syntax:
            %   maps = FittingProcessorRician_3Comp.FitImage(imData, sigma)
            %
            % Inputs:
            %   imData - structure with fields:
            %            .images         5-D array (nx,ny,nz,ncoils,nTE)
            %            .TE             1-by-n vector (echo times in seconds)
            %            .FieldStrength  Scalar (Tesla)
            %            .fittingIndent  (optional) indentation parameter
            %   sigma  - fixed noise sigma for fitting.
            %
            % Outputs:
            %   maps - structure with parameter maps including:
            %          FF (% fat), R2fat, R2water_long, R2water_short, T2* maps,
            %          S0, total Water, Fat, and individual water components.
            
            TE = 1000 * imData.TE;
            %TE = imData.TE;
            TE = reshape(TE, [], 1);
            FieldStrength = imData.FieldStrength;
            data = imData.images;
            [nx, ny, nz, nTE] = size(data);
            
            if isfield(imData, 'mask')
                mask = imData.mask;
            else
                Img2norm = data(:,:,:,1);
                threshold = prctile(Img2norm(:), 85);
                mask = Img2norm > threshold;
            end
            
            tmp = squeeze(sum(sum(mask,1),2));
            slice_idxs = find(tmp > 0);
            
            % Preallocate output vectors.
            FFrician = zeros(nx*ny*nz,1);
            R2fat_rician = zeros(nx*ny*nz,1);
            R2water_long_rician = zeros(nx*ny*nz,1);
            R2water_short_rician = zeros(nx*ny*nz,1);
            s0rician = zeros(nx*ny*nz,1);
            Wl_est = zeros(nx*ny*nz,1);
            Ws_est = zeros(nx*ny*nz,1);
            
            mask_r = mask(:);
            data_r = reshape(data, nx*ny*nz, nTE);
            RoiIdxs = mask_r > 0 & all(data_r > 0, 2);
            data_roi = data_r(RoiIdxs, :);
            nb_voxels = size(data_roi, 1);
            
            disp('Setting fit parameters for 3-component model...');
            tstart = tic;
            algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(data_roi, TE, sigma, 1);
            tstop = toc(tstart);
            disp("Parameters set in " + tstop + " seconds");
            
            % Preallocate results for each voxel.
            FFrician_roi = zeros(nb_voxels,1);
            R2fat_roi = zeros(nb_voxels,1);
            R2water_long_roi = zeros(nb_voxels,1);
            R2water_short_roi = zeros(nb_voxels,1);
            s0rician_roi = zeros(nb_voxels,1);
            Wl_roi = zeros(nb_voxels,1);
            Ws_roi = zeros(nb_voxels,1);
            
            disp('Fitting 3-component data...');
            % Set up a progress bar via DataQueue and waitbar.
            D = parallel.pool.DataQueue;
            progressCount = 0;
            hWaitbar = waitbar(0, 'Performing 3-component mapping...');
            afterEach(D, @updateWaitbar);
            function updateWaitbar(~)
                progressCount = progressCount + 1;
                waitbar(progressCount/nb_voxels, hWaitbar, ...
                    sprintf('Progress: %d%%', round(100*progressCount/nb_voxels)));
                drawnow;
            end
            
            tstart = tic;
            parfor i = 1:nb_voxels
                ydata = data_roi(i,:);
                algoparams = algoparamsArray(i);
                % Use the three-component fitting routine.
                outparams = FittingProcessorRician_3Comp.RicianMagnitudeFitting(TE, FieldStrength, ydata, algoparams.sigEst, algoparams);
                % The returned parameters: [F; Wl; Ws; R2fat; R2water_long; R2water_short]
                FFrician_roi(i) = outparams.F / (outparams.F + outparams.Wl + outparams.Ws);
                R2fat_roi(i) = outparams.R2fat;
                R2water_long_roi(i) = outparams.R2water_long;
                R2water_short_roi(i) = outparams.R2water_short;
                s0rician_roi(i) = outparams.F + outparams.Wl + outparams.Ws;
                Wl_roi(i) = outparams.Wl;
                Ws_roi(i) = outparams.Ws;
                send(D, i);
            end
            tstop = toc(tstart);
            delete(hWaitbar);
            disp("3-component data fitted in " + tstop + " seconds");
            
            % (No post-fit forcing is done; the water split remains as fitted.)
            
            FFrician(RoiIdxs,:) = FFrician_roi;
            R2fat_rician(RoiIdxs,:) = R2fat_roi;
            R2water_long_rician(RoiIdxs,:) = R2water_long_roi;
            R2water_short_rician(RoiIdxs,:) = R2water_short_roi;
            s0rician(RoiIdxs,:) = s0rician_roi;
            Wl_est(RoiIdxs,:) = Wl_roi;
            Ws_est(RoiIdxs,:) = Ws_roi;
            
            % Assemble output maps.
            maps.FF = reshape(100.*FFrician, nx, ny, nz);
            maps.R2fat = reshape(1000.*R2fat_rician, nx, ny, nz);
            maps.R2water_long = reshape(1000.*R2water_long_rician, nx, ny, nz);
            maps.R2water_short = reshape(1000.*R2water_short_rician, nx, ny, nz);
            maps.T2str_fat = 1./reshape(R2fat_rician, nx, ny, nz);
            maps.T2str_water_long = 1./reshape(R2water_long_rician, nx, ny, nz);
            maps.T2str_water_short = 1./reshape(R2water_short_rician, nx, ny, nz);
            maps.S0 = reshape(s0rician, nx, ny, nz);
            maps.Fat = reshape(s0rician .* FFrician, nx, ny, nz);
            maps.Water = reshape(s0rician - s0rician .* FFrician, nx, ny, nz);
            maps.Water_long = reshape(Wl_est, nx, ny, nz);
            maps.Water_short = reshape(Ws_est, nx, ny, nz);
            maps.Wsfraction = 100 .* maps.Water_short ./(maps.Fat + maps.Water_short + maps.Water_long);
            
            % Display a representative slice.
            if numel(slice_idxs) > 1
                slice = slice_idxs(round(numel(slice_idxs)/2));
            else
                slice = slice_idxs;
            end
            
            figure;
            subplot(2,5,1); imagesc(maps.FF(:,:,slice)); colorbar; title('FF (%)'); clim([0 100]);
            subplot(2,5,2); imagesc(maps.R2fat(:,:,slice)); colorbar; title('R2_{fat} (s^{-1})'); clim([0 250]);
            subplot(2,5,3); imagesc(maps.R2water_long(:,:,slice)); colorbar; title('R2_{water\_long} (s^{-1})'); clim([0 200]);
            subplot(2,5,4); imagesc(maps.R2water_short(:,:,slice)); colorbar; title('R2_{water\_short} (s^{-1})'); clim([0 200]);
            subplot(2,5,5); imagesc(maps.S0(:,:,slice)); colorbar; title('S0');
            subplot(2,5,6); imagesc(maps.T2str_fat(:,:,slice)); colorbar; title('T2*_{fat} (ms)'); clim([0 80]);
            subplot(2,5,7); imagesc(maps.T2str_water_long(:,:,slice)); colorbar; title('T2*_{water\_long} (ms)'); clim([0 80]);
            subplot(2,5,8); imagesc(maps.T2str_water_short(:,:,slice)); colorbar; title('T2*_{water\_short} (ms)'); clim([0 10]);
            subplot(2,5,9); imagesc(maps.Wsfraction(:,:,slice)); colorbar; title('Water\_short Fraction'); clim([0 100]);
        end
        
        function maps = MultistepFitImage(imData, roi)
            % MultistepFitImage combines sigma estimation and fixed-sigma fitting.
            % It first estimates sigma using the 2-component model and then performs the
            % 3-component fit.
            %
            % Syntax:
            %   maps = FittingProcessorRician_3Comp.MultistepFitImage(imData, roi)
            %
            % Inputs:
            %   imData - structure with fields: .images, .TE, .FieldStrength, .fittingIndent
            %   roi    - structure with field: .mask (nx-by-ny logical)
            %
            % Outputs:
            %   maps - Final parameter maps.
            
            [~, ~, sigmaFromFit] = FittingProcessorRician_3Comp.FitImageSigma(imData, roi);
            disp("Estimated sigma is: " + sigmaFromFit);
            mapsFittedSigma = FittingProcessorRician_3Comp.FitImage(imData, sigmaFromFit);
            maps = mapsFittedSigma;
        end
        
        function algoparamsArray = setAlgoparamsArray(S, echotimes, sigmaEstimateFromRoi, opt)
            % setAlgoparamsArray creates an array of fitting parameter structures for each voxel.
            %
            % For the 3-component model (opt == 1), the parameter vector is:
            %   [F; Wl; Ws; R2fat; R2water_long; R2water_short]
            %
            % For sigma estimation (opt == 2), the standard 2-component model is used:
            %   [F; W; R2water; R2fat; sigma]
            %
            % Syntax:
            %   algoparamsArray = FittingProcessorRician_3Comp.setAlgoparamsArray(S, echotimes, sigmaEstimateFromRoi, opt)
            
            nvoxels = size(S, 1);
            solver = 'fmincon';
            opts = optimoptions('fmincon', ...
                'Algorithm', 'interior-point', ...
                'InitBarrierParam', 100000, ...
                'ScaleProblem', true, ...
                'FiniteDifferenceType', 'central', ...
                'Display', 'off');
            
            if opt == 1
                % 3-component model.
                algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                    'lb', [], 'ub', [], 'solver', solver, 'options', opts), nvoxels, 1);
            elseif opt == 2
                % 2-component model for sigma estimation.
                algoparamsArray = repmat(struct('vinit', [], 'Sinit', [], 'sigEst', [], ...
                    'lb', [], 'ub', [], 'solver', solver, 'options', opts), nvoxels, 1);
            else
                error('Unknown option value: %d. Valid options are 1 or 2', opt);
            end
            
            vinit = 0.1;    % initial guess for R2 (for R2fat and R2water_long)
            vmax = 5;       % maximum R2 value (ms^-1)
            vmin = 0.5e-3;  % minimum R2 value (ms^-1)
            [maxS, idx] = max(S, [], 2);
            Sinit = maxS .* exp(vinit * echotimes(idx));
            
            if opt == 1
                % For 3-component model, split water equally initially.
                parfor i = 1:nvoxels
                    % Use half-and-half for water initialization.
                    % Parameter vector: [F; Wl; Ws; R2fat; R2water_long; R2water_short]
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    % % Lower bound for R2water_short is set to 1.
                    % algoparamsArray(i).lb = [0; 0; 0; vmin; vmin; 1];
                    % algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); 3 * Sinit(i); vmax; vmax; vmax];
                    % Lower bound for R2wshort is 0.167, limiting T2*w,s at
                    % 6 ms; R2wlong bound from 6.01-20 ms
                    %algoparamsArray(i).lb = [0; 0; 0; vmin; 0.05; 0.166];
                    %algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); 3 * Sinit(i); vmax; 0.167; 100];
                    % limit at 4 ms; R2wlong bound from 4.01-20 ms
                    algoparamsArray(i).lb = [0; 0; 0; vmin; 0.05; 0.249];
                    algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); 3 * Sinit(i); vmax; 0.25; 100];
                    % lower bound for R2wshort is 0.125 (T2*ws at 8 ms)
                    %algoparamsArray(i).lb = [0; 0; 0; vmin; 0.05; 0.1249];
                    %algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); 3 * Sinit(i); vmax; 0.125; 100];
                    % same limits see what happens
                    %algoparamsArray(i).lb = [0; 0; 0; vmin; vmin; vmin];
                    %algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); 3 * Sinit(i); vmax; vmax; 100];
                end
            elseif opt == 2
                % 2-component model.
                parfor i = 1:nvoxels
                    algoparamsArray(i).vinit = vinit;
                    algoparamsArray(i).Sinit = Sinit(i);
                    algoparamsArray(i).sigEst = sigmaEstimateFromRoi;
                    algoparamsArray(i).lb = [0; 0; vmin; vmin; sigmaEstimateFromRoi*0.2];
                    algoparamsArray(i).ub = [3 * Sinit(i); 3 * Sinit(i); vmax; vmax; sigmaEstimateFromRoi*2];
                end
            end
        end
        
        function results = RicianMagnitudeFitting(echotimes, tesla, Smagnitude, sig, algoparams)
            % RicianMagnitudeFitting performs Rician magnitude fitting with multiple
            % initializations for the 3-component model.
            %
            % Parameter vector: [F; Wl; Ws; R2fat; R2water_long; R2water_short]
            %
            % Syntax:
            %   results = FittingProcessorRician_3Comp.RicianMagnitudeFitting(echotimes, tesla, Smagnitude, sig, algoparams)
            
            RicianFit.solver = algoparams.solver;
            RicianFit.options = algoparams.options;
            RicianFit.lb = algoparams.lb;
            RicianFit.ub = algoparams.ub;
            RicianFit.objective = @(p) -FittingProcessorRician_3Comp.R2RicianObj(p, echotimes, tesla, Smagnitude, sig);
            
            % Define several initial guesses that enforce half-and-half splitting in water.
            % modified the initial guesses from MB version
            initialGuesses = { [0.001; algoparams.Sinit/2; algoparams.Sinit/2; algoparams.vinit; algoparams.vinit; 0.2], ...
                               [algoparams.Sinit; 0.001; 0.001; algoparams.vinit; algoparams.vinit; 0.2]};
            nInit = numel(initialGuesses);
            pminCell = cell(nInit,1);
            fminArray = zeros(nInit,1);
            for k = 1:nInit
                RicianFit.x0 = initialGuesses{k};
                [p, fval] = fmincon(RicianFit);
                pminCell{k} = p;
                fminArray(k) = fval;
            end
            [results.fmin, idx] = min(fminArray);
            results.chosenmin = idx;
            bestP = pminCell{results.chosenmin};
            results.F = bestP(1);
            results.Wl = bestP(2);
            results.Ws = bestP(3);
            results.R2fat = bestP(4);
            results.R2water_long = bestP(5);
            results.R2water_short = bestP(6);
            results.SSE = FittingProcessorRician_3Comp.GetSEE(bestP, echotimes, tesla, Smagnitude);
        end
        
        function results = RicianMagnitudeFitting_WithSigma(echotimes, tesla, Smagnitude, sigmaInit, algoparams)
            % RicianMagnitudeFitting_WithSigma performs Rician fitting with sigma estimation.
            % For sigma estimation the 2-component model is used.
            %
            % Parameter vector: [F; W; R2water; R2fat; sigma]
            %
            % Syntax:
            %   results = FittingProcessorRician_3Comp.RicianMagnitudeFitting_WithSigma(echotimes, tesla, Smagnitude, sigmaInit, algoparams)
            
            Ricianfitting.solver = algoparams.solver;
            Ricianfitting.options = algoparams.options;
            Ricianfitting.lb = algoparams.lb;
            Ricianfitting.ub = algoparams.ub;
            Ricianfitting.objective = @(p) -FittingProcessorRician_3Comp.R2RicianObj_WithSigma(p, echotimes, tesla, Smagnitude);
            
            Ricianfitting.x0 = [0.001; algoparams.Sinit; algoparams.vinit; algoparams.vinit; sigmaInit];
            [pmin1, fmin1] = fmincon(Ricianfitting);
            Ricianfitting.x0 = [algoparams.Sinit; 0.001; algoparams.vinit; algoparams.vinit; sigmaInit];
            [pmin2, fmin2] = fmincon(Ricianfitting);
            
            if fmin1 <= fmin2
                results.F = pmin1(1);
                results.W = pmin1(2);
                results.R2water = pmin1(3);
                results.R2fat = pmin1(4);
                results.sig = pmin1(5);
                results.fmin = fmin1;
                results.chosenmin = 1;
                results.SSE = FittingProcessorRician_3Comp.GetSEE(pmin1, echotimes, tesla, Smagnitude);
            else
                results.F = pmin2(1);
                results.W = pmin2(2);
                results.R2water = pmin2(3);
                results.R2fat = pmin2(4);
                results.sig = pmin2(5);
                results.fmin = fmin2;
                results.chosenmin = 2;
                results.SSE = FittingProcessorRician_3Comp.GetSEE(pmin2, echotimes, tesla, Smagnitude);
            end
        end
        
        function sse = GetSEE(p, echotimes, tesla, Smeasured)
            % If using the 3-component model, p should have 6 elements.
            % If using the 2-component model (for sigma estimation), p will have 5.
            if numel(p) == 6
                Spredicted = FittingProcessorRician_3Comp.MultiPeakFatTripleR2(echotimes, tesla, ...
                    p(1), p(2), p(3), p(4), p(5), p(6), 0);
            elseif numel(p) == 5
                % Use the two-component model: parameters [F; W; R2water; R2fat; sigma]
                Spredicted = FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(echotimes, tesla, ...
                    p(1), p(2), p(3), p(4), 0);
            else
                error('Parameter vector has insufficient elements.');
            end
            errors = abs(Spredicted) - abs(Smeasured);
            sse = sum(errors.^2);
        end

        
        function [loglik] = R2RicianObj(p, echotimes, tesla, Smeasured, sig)
            % R2RicianObj computes the log likelihood for Rician fitting for the 3-component model.
            %
            % Syntax:
            %   loglik = FittingProcessorRician_3Comp.R2RicianObj(p, echotimes, tesla, Smeasured, sig)
            
            Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatTripleR2(echotimes, tesla, p(1), p(2), p(3), p(4), p(5), p(6), 0));
            Smeasured = abs(Smeasured);
            loglik = FittingProcessorRician_3Comp.RicianLogLik(Smeasured, Spredicted, sig);
        end
        
        function [loglik] = R2RicianObj_WithSigma(p, echotimes, tesla, Smeasured)
            % R2RicianObj_WithSigma computes the log likelihood for Rician fitting with sigma estimation.
            % This function uses the 2-component model.
            %
            % Syntax:
            %   loglik = FittingProcessorRician_3Comp.R2RicianObj_WithSigma(p, echotimes, tesla, Smeasured)
            
            Spredicted = abs(FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(echotimes, tesla, p(1), p(2), p(3), p(4), 0));
            Smeasured = abs(Smeasured);
            loglik = FittingProcessorRician_3Comp.RicianLogLik(Smeasured, Spredicted, p(5));
        end
        
        function S = MultiPeakFatTripleR2(t, tesla, F, Wl, Ws, R2fat, R2water_long, R2water_short, fB)
            % MultiPeakFatTripleR2 computes the complex MRI signal from fat and water
            % using a three-component model. The signal is computed as the sum of:
            %   fatComponent = F * fatSignal * exp(-R2fat*t)
            %   waterLongComponent = Wl * waterSignal * exp(-R2water_long*t)
            %   waterShortComponent = Ws * waterSignal * exp(-R2water_short*t)
            % Off-resonance is then applied.
            %
            % Syntax:
            %   S = FittingProcessorRician_3Comp.MultiPeakFatTripleR2(t, tesla, F, Wl, Ws, R2fat, R2water_long, R2water_short, fB)
            %
            % Inputs:
            %   t             - Vector of echo times in milliseconds.
            %   tesla         - Field strength in Tesla.
            %   F             - Fat signal at time zero.
            %   Wl            - Water_long signal at time zero.
            %   Ws            - Water_short signal at time zero.
            %   R2fat         - Decay constant for fat.
            %   R2water_long  - Decay constant for water_long.
            %   R2water_short - Decay constant for water_short.
            %   fB            - Off-resonance frequency in Hz.
            %
            % Output:
            %   S             - Complex-valued signal.
            
            gyro = 42.58e6;
            larmor = tesla * gyro;
            
            % Fat frequency parameters.
            Fatamps = [0.087; 0.694; 0.128; 0.004; 0.039; 0.048];
            Fatshift = [-3.9; -3.5; -2.7; -2.04; -0.49; 0.50];
            
            % Convert frequency shifts from ppm to rad/ms.
            FatFreqs = Fatshift * larmor * 1e-6;
            FatFreqs = FatFreqs / 1000;
            Fatw = 2*pi*FatFreqs;
            
            Waterw = 0;
            
            fatExp = exp(1i * Fatw * t');  % 6-by-T
            fatSignal = Fatamps.' * fatExp;  % 1-by-T
            
            waterExp = exp(1i * Waterw * t');  % 1-by-T
            
            fatComponent = F * fatSignal .* exp(-R2fat*t');
            waterLongComponent = Wl * waterExp .* exp(-R2water_long*t');
            waterShortComponent = Ws * waterExp .* exp(-R2water_short*t');
            
            S = fatComponent + waterLongComponent + waterShortComponent;
            
            fB = fB/1000;
            offResonance = exp(1i * 2*pi*fB*t');
            S = S .* offResonance;
        end
        
        function S = MultiPeakFatDoubleR2(t, tesla, F, W, R2water, R2fat, fB)
            % MultiPeakFatDoubleR2 computes the complex MRI signal from water and fat
            % using a two-component model. This function is used for sigma estimation.
            %
            % Syntax:
            %   S = FittingProcessorRician_3Comp.MultiPeakFatDoubleR2(t, tesla, F, W, R2water, R2fat, fB)
            %
            % Inputs:
            %   t       - Vector of echo times in milliseconds.
            %   tesla   - Field strength in Tesla.
            %   F       - Fat signal at time zero.
            %   W       - Water signal at time zero.
            %   R2water - Decay constant for water.
            %   R2fat   - Decay constant for fat.
            %   fB      - Off-resonance frequency in Hz.
            %
            % Output:
            %   S       - Complex-valued signal.
            
            gyro = 42.58e6;
            larmor = tesla * gyro;
            
            Fatamps = [0.087; 0.694; 0.128; 0.004; 0.039; 0.048];
            Fatshift = [-3.9; -3.5; -2.7; -2.04; -0.49; 0.50];
            
            FatFreqs = Fatshift * larmor * 1e-6;
            FatFreqs = FatFreqs/1000;
            Fatw = 2*pi*FatFreqs;
            
            Waterw = 0;
            
            fatExp = exp(1i * Fatw * t');
            fatSignal = Fatamps.' * fatExp;
            
            waterExp = exp(1i * Waterw * t');
            
            fatComponent = F * fatSignal .* exp(-R2fat*t');
            waterComponent = W * waterExp .* exp(-R2water*t');
            
            S = fatComponent + waterComponent;
            
            fB = fB/1000;
            offResonance = exp(1i * 2*pi*fB*t');
            S = S .* offResonance;
        end
        
        function [loglik, logliks] = RicianLogLik(measurements, predictions, sigma)
            % RicianLogLik computes the log likelihood for the Rician noise model.
            sigmaSquared = sigma.^2;
            sumsqsc = (measurements.^2 + predictions.^2) ./ (2*sigmaSquared);
            scp = measurements.*predictions ./ sigmaSquared;
            lb0 = FittingProcessorRician_3Comp.logbesseli0(scp);
            logliks = log(measurements) - log(sigmaSquared) - sumsqsc + lb0;
            loglik = sum(logliks);
        end
        
        function lb0 = logbesseli0(x)
            % logbesseli0 computes a robust logarithm of besseli(0,x).
            lb0 = zeros(size(x));
            exact = x < 700;
            approx = ~exact;
            lb0(exact) = log(besseli(0,x(exact)));
            lb0(approx) = x(approx) - 0.5*log(2*pi*x(approx));
        end
        
        function loglik = GaussianLogLik(measured, predicted, sigma)
            % GaussianLogLik computes the Gaussian log likelihood.
            errors = measured - predicted;
            n = numel(errors);
            loglik = -0.5*sum((errors./sigma).^2) - n*log(sigma*sqrt(2*pi));
        end
        
        function img_whiten = WhitenImage(img)
            % WhitenImage performs quartile normalization on a 3D image.
            img_whiten = img;
            num_slices = size(img,3);
            for s = 1:num_slices
                slice = img(:,:,s);
                centered = slice - median(slice(:));
                normalized = centered / prctile(centered(:),75);
                p25 = prctile(normalized(:),25);
                p75 = prctile(normalized(:),75);
                scaleFactor = 2 / (p75 - p25);
                offset = -1 - 2*p25/(p75-p25);
                img_clean = normalized*scaleFactor + offset;
                lower_bound = prctile(img_clean(:),3);
                upper_bound = prctile(img_clean(:),97);
                img_clean = max(min(img_clean,upper_bound),lower_bound);
                img_whiten(:,:,s) = img_clean;
            end
        end
    end
end
