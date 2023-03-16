% Clear variables
close all; clearvars; clc;

outputImages = true;

% Create paths to load and store data from
path = pwd;
dpath = strrep(path, 'data', 'code/experiments');

% Ensure data folders exists
if ~exist([dpath '/Detection validation'], 'dir')
    mkdir([dpath '/Detection validation'])
end

if ~exist([dpath '/Processed data'], 'dir')
    mkdir([dpath '/Processed data'])
end

if ~exist([dpath '/Radius Curves'], 'dir')
    mkdir([dpath '/Radius Curves'])
end

if ~exist([dpath '/Example'], 'dir')
    mkdir([dpath '/Example'])
end

% number of colonies and time steps
nT = 33;
nC = 30;
dp = 0.65; % Size of pixel (in micron)

% Time range
TT = [5 (0:nT-2)*0.5+6+1/6];
B = [1 1 1 1 1 3 3 3 3 3 5 5 5 5 5 2 2 2 2 2 4 4 4 4 4 6 6 6 6 6];
dT = (5-min(B, 5))*0.5;
T = TT' + dT;

% Detmine labels
XX = repmat(arrayfun(@(k)sum(B(1:k)==B(k)), 1:nC), nT, 1);
BB = repmat(B, nT, 1);

% First determine the BF radius of all images %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% while doing so, determine the GFP expression inside the BF radius and
% outside

% Prepare radius array
radius  = nan(nT, nC);

% Determine the BG and max level of GFP images
BG_GFP  = zeros(nT, nC);
MAX_GFP = zeros(nT, nC);


% Loop over series
for k = 1:nT*nC

    s = ceil(k/ nT);
    t = mod(k-1, nT) + 1;

    % Check the data exists
    if ~exist([path '/Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c001.tif', s, t)], 'file')
        continue
    end

    % Load the picture
    data  = double(importdata([path '/Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c001.tif', s, t)]));
    gdata = double(importdata([path '/Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c002.tif', s, t)]));

    % Use median filter to smooth data
    fdata = medfilt2(data);

    % Remove border
    m = median(fdata(:));
    [X, Y] = meshgrid(1:size(fdata, 1), 1:size(fdata, 2));
    d = sqrt((X - 1150).^2 + (Y - 1330).^2);
    fdata(d'>1250) = m;

    % Store the min/max values
    level_black_BF = 100*min(fdata(:))/2^16;
    level_white_BF = 100*max(fdata(:))/2^16;

    % Create structue element
    se = strel('disk', 50);

    fdata = imopen(fdata, se);
    fdata = imclose(fdata, se);

    % Rescale data
    fdata = 1 - fdata/2^16;
    fdata = (fdata-median(fdata(:)))/(1-median(fdata(:)));

    % Do some more filtering on the BF data
    xx = 1:size(fdata, 2);
    yy = 1:size(fdata, 1);
    ff = fdata;

    % Adjust for background gradient
    ff(d'>1250) = nan;
    ff(d'<1100) = nan;

    [x, y] = meshgrid(xx, yy);
    f = fit([x(:), y(:)], fdata(:), 'poly11');
    fdata = fdata-f(x, y);

    % Apply the theshold
    bw = fdata>0.2;
    bw = imclose(bw, se);

    % Determine the largest region
    L = bwlabel(bw);
    areas = arrayfun(@(d)bwarea(L==d), 1:max(L(:)));
    [~, I] = sort(areas, 'desc');

    % Filter out regions touching the edge
    for i = I
        bw = L==i;

        % Check if the edges are hit
        if ~any([bw(1, :) bw(end, :) bw(:, 1)' bw(:, end)'])
            break
        end

        % If nothing is found, zero the bw
        if i == I(end)
            bw = zeros(size(bw));
        end
    end

    % Get image properties of largest cells
    stat  = regionprops(bw, 'Area', 'EquivDiameter', 'Centroid');

    % Store radius and mean signals
    if ~isempty(stat)
        radius(k)  = stat.EquivDiameter/2*dp;
        BG_GFP(k)  = median(gdata(~bw(:)));
        MAX_GFP(k) = median(gdata(bw(:)));
    end

    if outputImages

        % Detect boundary
        bfBorder = zeros(size(bw));
        for i = 1:size(bw, 1)
            for j = 1:size(bw, 2)
                if bw(i, j) == 1
                    bfBorder(i, j) = 1;
                    continue;
                end
                for ii = max(1, i-1):min(size(bw, 1), i+1)
                    for jj = max(1, j-1):min(size(bw, 2), j+1)
                        if bw(ii,jj) == 1
                            bfBorder(i, j) = 1;
                        end
                    end
                end
            end
        end
        bfBorder = bfBorder - bw;


        % Process the BF images
        inputPath = ['Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c001.tif', s, t)];
        tempPath  = ['Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c001_temp.tif', s, t)];
        outputPath1 = ['../code/experiments/Processed data/' sprintf('B%dX%d_t%0.3d_c001.tif', BB(k), XX(k), t)];
        outputPath2 = ['../code/experiments/Detection validation/' sprintf('B%dX%d_t%0.3d_c001.tif', BB(k), XX(k), t)];

        % Fix the leveling of the BF image
        cmd = sprintf('convert "%s" -level %.3f%%,%.3f%% -depth 8 -type TrueColor "%s"', inputPath, level_black_BF, level_white_BF, tempPath);
        system([cmd ' 2> /dev/null']);

        % Overlay the detection lines
        cmd = sprintf('cp "%s" "%s"', tempPath, outputPath2);
        system(cmd);

        cmdStart = sprintf('convert "%s" -fill "rgb(255,0,0)" ', outputPath2);
        cmd = [];
        for i = 1:size(bfBorder, 1)
            for j = 1:size(bfBorder, 2)
                if bfBorder(i,j) == 1

                    if isempty(cmd)
                        cmd = cmdStart;
                    end

                    cmd = [cmd sprintf('-draw "color %d,%d point" ', j, i)];
                    if numel(cmd) > 500
                        cmd = [cmd sprintf('"%s"', outputPath2)];
                        err = system(cmd);
                        if err == 1
                            cmd
                        end
                        cmd = [];
                    end
                end
            end
        end
        if ~isempty(cmd)
            cmd = [cmd sprintf('"%s"', outputPath2)];
            err = system(cmd);
            if err == 1
                cmd
            end
        end

        % Resize the output of BF images
        cmd = sprintf('convert "%s" -resize 25%% "%s"', outputPath2, outputPath2);
        system(cmd);

        % Annotate the clean BF images
        cmd = sprintf('convert "%s" -resize 25%% "%s"', tempPath, outputPath1);
        system(cmd);

        cmd = sprintf('convert "%s" -gravity north -stroke none -fill white -pointsize 50 -annotate 0 "%s" "%s"', outputPath1, sprintf('BF, T = %.1f h', T(k)), outputPath1);
        system(cmd);

        cmd = sprintf('convert "%s" -gravity southeast -stroke none -fill white -pointsize 36 -annotate +5+20 "100 µm" "%s"', outputPath1, outputPath1);
        system(cmd);

        cmd = sprintf('composite -compose over -gravity southeast -geometry +5+5 \\( -size %dx10 xc:white \\) "%s" "%s"', round(100/dp), outputPath1, outputPath1);
        system(cmd);

        % Clean up
        cmd = sprintf('rm "%s"', tempPath);
        system(cmd);

        % Create special output for manuscript figure
        if (any([and(BB(k) == 5, XX(k) == 4) and(BB(k) == 3, XX(k) == 5) and(BB(k) == 2, XX(k) == 4) and(BB(k) == 1, XX(k) == 1) and(BB(k) == 1, XX(k) == 5) and(BB(k) == 3, XX(k) == 2)]) && any(t == 4:2:24))

            % Process the BF images
            inputPath  = ['Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c001.tif', s, t)];
            outputPath = ['../code/experiments/Example/' sprintf('B%dX%d_t%0.3d_BF.tif', BB(k), XX(k), t)];

            % Fix the leveling of the BF image
            cmd = sprintf('convert "%s" -level %.3f%%,%.3f%% -depth 8 -type TrueColor "%s"', inputPath, level_black_BF, level_white_BF, outputPath);
            system([cmd ' 2> /dev/null']);

        end
    end

end


% Next determine the GFP radius of all images %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Prepare radius array
thresholds = unique([1:0.1:2.5 1.8:0.025:2.0]); % Adding extra fidelity around where \nu seems to be

GFP_radius = nan(nT, nC);
for threshold = thresholds
    eval(['GFP_radius_' strrep(sprintf('%.2f', threshold), '.', '_') ' = nan(nT, nC)'])
end

% Use control colonies for calibration
maxGFP = median(MAX_GFP(:, end-4:end), 2);
bgGFP  = median(BG_GFP(:, end-4:end), 2);

% Parfor recasting (prevent broadcast variable?)
maxGFP = repmat(maxGFP, 1, nC);
bgGFP  = repmat(bgGFP, 1, nC);

% Loop over series
for k = 1:nT*nC

    s = ceil(k/ nT);
    t = mod(k-1, nT) + 1;

    %     medianGFP_t = bgGFP(k);
    %     maxGFP_t    = maxGFP(k);

    medianGFP_t = BG_GFP(k);
    maxGFP_t    = maxGFP(k);

    % Check the data exists
    if ~exist([path '/Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c002.tif', s, t)], 'file')
        continue
    end

    % Load the picture
    gdata = double(importdata([path '/Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c002.tif', s, t)]));

    % Create structue element
    se = strel('disk', 50);

    % Detect the GFP radius
    fdata = imopen(gdata, se);
    fdata = imclose(fdata, se);
    
    % Loop over thresholds
    for threshold = thresholds

        % Apply the threshold
        bw = fdata > threshold * medianGFP_t;
        bw = imclose(bw, se);

        % Determine the largest region
        L = bwlabel(bw);
        areas = arrayfun(@(d)bwarea(L==d), 1:max(L(:)));
        [~, I] = sort(areas, 'desc');

        % Filter out regions touching the edge
        for i = I
            bw = L==i;

            % Check if the edges are hit
            if ~any([bw(1, :) bw(end, :) bw(:, 1)' bw(:, end)'])
                break
            end

            % If nothing is found, zero the bw
            if i == I(end)
                bw = zeros(size(bw));
            end
        end

        % Get image properties of largest cells
        stat  = regionprops(bw, 'Area', 'EquivDiameter', 'Centroid');

        % Store radius and mean signals
        if ~isempty(stat)
            eval(['GFP_radius_' strrep(sprintf('%.2f', threshold), '.', '_') '(k) = stat.EquivDiameter/2*dp;'])
        end

        if outputImages && round(threshold, 2) == 1.95
            % Detect boundary
            gfpBorder = zeros(size(bw));
            for i = 1:size(bw, 1)
                for j = 1:size(bw, 2)
                    if bw(i, j) == 1
                        gfpBorder(i, j) = 1;
                        continue;
                    end
                    for ii = max(1, i-1):min(size(bw, 1), i+1)
                        for jj = max(1, j-1):min(size(bw, 2), j+1)
                            if bw(ii,jj) == 1
                                gfpBorder(i, j) = 1;
                            end
                        end
                    end
                end
            end
            gfpBorder = gfpBorder - bw;
        end
    end

    if outputImages

        % Process the GFP images
        inputPath = ['Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c002.tif', s, t)];
        tempPath  = ['Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c002_temp.tif', s, t)];
        outputPath1 = ['../code/experiments/Processed data/' sprintf('B%dX%d_t%0.3d_c002.tif', BB(k), XX(k), t)];
        outputPath2 = ['../code/experiments/Detection validation/' sprintf('B%dX%d_t%0.3d_c002.tif', BB(k), XX(k), t)];

        level_black_GFP = 100*medianGFP_t/2^16;
        level_white_GFP = 100*5*medianGFP_t/2^16;

        % Fix the leveling of the GFP image
        cmd = sprintf('convert "%s" -level %.3f%%,%.3f%% -fill green -tint 100 "%s"', inputPath, level_black_GFP, level_white_GFP, tempPath);
        system([cmd ' 2> /dev/null']);

        % Overlay the detection lines
        cmd = sprintf('cp "%s" "%s"', tempPath, outputPath2);
        system(cmd);

        cmdStart = sprintf('convert "%s" -fill "rgb(255,0,0)" ', outputPath2);
        cmd = [];
        for i = 1:size(gfpBorder, 1)
            for j = 1:size(gfpBorder, 2)
                if gfpBorder(i,j) == 1

                    if isempty(cmd)
                        cmd = cmdStart;
                    end

                    cmd = [cmd sprintf('-draw "color %d,%d point" ', j, i)];
                    if numel(cmd) > 500
                        cmd = [cmd sprintf('"%s"', outputPath2)];
                        err = system(cmd);
                        if err == 1
                            cmd
                        end
                        cmd = [];
                    end
                end
            end
        end
        if ~isempty(cmd)
            cmd = [cmd sprintf('"%s"', outputPath2)];
            err = system(cmd);
            if err == 1
                cmd
            end
        end

        % Resize the output of GFP images
        cmd = sprintf('convert "%s" -resize 25%% "%s"', outputPath2, outputPath2);
        system(cmd);


        % Annotate the clean GFP images
        cmd = sprintf('convert "%s" -resize 25%% "%s"', tempPath, outputPath1);
        system(cmd);

        cmd = sprintf('convert "%s" -gravity north -stroke none -fill white -pointsize 50 -annotate 0 "%s" "%s"', outputPath1, sprintf('GFP, T = %.1f h', T(k)), outputPath1);
        system(cmd);

        cmd = sprintf('convert "%s" -gravity southeast -stroke none -fill white -pointsize 36 -annotate +5+20 "100 µm" "%s"', outputPath1, outputPath1);
        system(cmd);

        cmd = sprintf('composite -compose over -gravity southeast -geometry +5+5 \\( -size %dx10 xc:white \\) "%s" "%s"', round(100/dp), outputPath1, outputPath1);
        system(cmd);

        % Combine into comparrison image
        inputPath =@(c) ['../code/experiments/Processed data/' sprintf('B%dX%d_t%0.3d_c00%d.tif', BB(k), XX(k), t, c)];
        outputPath = ['../code/experiments/Processed data/' sprintf('B%dX%d_t%0.3d.tif', BB(k), XX(k), t)];

        cmd = sprintf('convert +append "%s" "%s" -resize 50%% "%s"', inputPath(1), inputPath(2), outputPath);
        system(strrep(cmd, '%', '%%'));

        % Clean up
        cmd = sprintf('rm "%s" "%s"', inputPath(1), inputPath(2));
        system(cmd);

        cmd = sprintf('rm "%s"', tempPath);
        system(cmd);

        % Create special output for manuscript figure
        if (any([and(BB(k) == 1, XX(k) == 1) and(BB(k) == 1, XX(k) == 4) and(BB(k) == 2, XX(k) == 2) and(BB(k) == 3, XX(k) == 5) and(BB(k) == 4, XX(k) == 5)]) && any(t == 4:2:24))

            % Process the GFP images
            inputPath = ['Exported data/Data_Post_Phage' sprintf('_Series_%d_t%0.3d_c002.tif', s, t)];
            outputPath = ['../code/experiments/Example/' sprintf('B%dX%d_t%0.3d_GFP.tif', BB(k), XX(k), t)];

            % Fix the leveling of the GFP image
            cmd = sprintf('convert "%s" -level %.3f%%,%.3f%% -fill green -tint 100 "%s"', inputPath, level_black_GFP, level_white_GFP, outputPath);
            system([cmd ' 2> /dev/null']);

        end
    end
end

save([dpath '/data.mat'], 'TT', 'dT', 'BG_GFP', 'MAX_GFP', 'radius', 'GFP_radius', 'GFP_radius_1_00', 'GFP_radius_1_10', 'GFP_radius_1_20', 'GFP_radius_1_30', 'GFP_radius_1_40', 'GFP_radius_1_50', 'GFP_radius_1_60', 'GFP_radius_1_70', 'GFP_radius_1_80', 'GFP_radius_1_90', 'GFP_radius_2_00', 'GFP_radius_2_10', 'GFP_radius_2_20', 'GFP_radius_2_30', 'GFP_radius_2_40', 'GFP_radius_2_50')
