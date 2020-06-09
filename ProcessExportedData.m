% Clear variables
close all; clearvars; clc;

showProgress = false;

% Create paths to load and store data from
path = pwd;
dpath = strrep(path, 'data', 'code/experiments');

% Ensure data folders exists
if ~exist([dpath '/Processed data'], 'dir')
    mkdir([dpath '/Processed data'])
end

if ~exist([dpath '/Radius Curves'], 'dir')
    mkdir([dpath '/Radius Curves'])
end

% Create output file
fid = fopen([path '/ConvertImages.sh'], 'w');

% number of colonies and time steps
nT = 33;
nC = 30;
dp = 0.65; % Size of pixel (in micron)

% Time range
TT = [5 (0:nT-1)*0.5+6+1/6];
C = [1 1 1 1 1 3 3 3 3 3 5 5 5 5 5 2 2 2 2 2 4 4 4 4 4 6 6 6 6 6];
dT = (5-min(C, 5))*0.5;


% Prepare radius array
radius = zeros(nT, nC);
GFP    = zeros(nT, nC);
GFP_radius_0025 = zeros(nT, nC);
GFP_radius_0030 = zeros(nT, nC);
GFP_radius_0035 = zeros(nT, nC);
GFP_radius_0040 = zeros(nT, nC);

% Determine the BG and max level of GFP images
BG_GFP = zeros(nT, nC);
MAX_GFP = zeros(nT, nC);
for s = 1:nC
    for k = 1:2
        switch k
            case 1
                T = 1;
                prefix = 'Data_Pre_Phage';
            case 2
                T = 2:34; %36
                prefix = 'Data_Post_Phage';
        end
        for t = 1:numel(T)
            A = imread([path '/Exported data/' sprintf('%s_Series_%d_t%0.3d_c002.tif', prefix, s, t)]);
            BG_GFP(t, s) = median(A(:));
            MAX_GFP(t, s) = max(A(:));
        end
    end
end

% Loop over series
for s = 1:nC

    % Loop over data series
    for k = 1:2
        switch k
            case 1
                T = 1;
                prefix = 'Data_Pre_Phage';
            case 2
                T = 2:34; %36
                prefix = 'Data_Post_Phage';
        end

        % Loop over time points before phage
        for t = 1:numel(T)

            % Check the data exists
            if ~exist([path '/Exported data/' sprintf('%s_Series_%d_t%0.3d_c001.tif', prefix, s, t)], 'file')
                continue
            end

            % Load the picture
            data  = double(importdata([path '/Exported data/' sprintf('%s_Series_%d_t%0.3d_c001.tif', prefix, s, t)]));
            gdata = double(importdata([path '/Exported data/' sprintf('%s_Series_%d_t%0.3d_c002.tif', prefix, s, t)]));

            if showProgress
                figure(1)
                imagesc(data)
                drawnow;
            end

            % Use median filter to smooth data
            fdata = medfilt2(data);

            % Remove border
            m = median(fdata(:));
            [X, Y] = meshgrid(1:size(fdata, 1), 1:size(fdata, 2));
            d = sqrt((X - 1150).^2 + (Y - 1330).^2);
            fdata(d'>1250) = m;

            % Store the min/max values
            level_black = 100*min(fdata(:))/2^16;
            level_white = 100*max(fdata(:))/2^16;

            % Create structue element
            se = strel('disk', 50);

            fdata = imopen(fdata, se);
            fdata = imclose(fdata, se);

            % Rescale data
            fdata = 1 - fdata/2^16;
            fdata = (fdata-median(fdata(:)))/(1-median(fdata(:)));

            if showProgress
                figure(1)
                imagesc(fdata)
                drawnow;
            end

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

            % Remove edge boundaries
            if showProgress
                figure(2)
                imagesc(bw)
                drawnow;
            end

            % Get image properties of largest cells
            stat  = regionprops(bw, 'Area', 'EquivDiameter', 'Centroid');

            % Store radius and mean signals
            if ~isempty(stat)
                radius(T(t), s) = stat.EquivDiameter/2*dp;
                GFP(T(t), s) = ((sum(sum(gdata(bw==1)))/sum(bw(:)))-median(BG_GFP(:)))/(max(MAX_GFP(:))-median(BG_GFP(:)));
            end

            % Detect the GFP radius
            fdata = imopen(gdata, se);
            fdata = imclose(fdata, se);

            % Rescale data
            fdata = fdata/2^16;

             % Apply the theshold
            bw = fdata>0.0025;
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
                GFP_radius_0025(T(t), s) = stat.EquivDiameter/2*dp;
            end

             % Apply the theshold
            bw = fdata>0.0030;
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
                GFP_radius_0030(T(t), s) = stat.EquivDiameter/2*dp;
            end

             % Apply the theshold
            bw = fdata>0.0035;
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
                GFP_radius_0035(T(t), s) = stat.EquivDiameter/2*dp;
            end

             % Apply the theshold
            bw = fdata>0.0040;
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
                GFP_radius_0040(T(t), s) = stat.EquivDiameter/2*dp;
            end

            % Fix the leveling of the BF image
            inputPath = ['Exported data/' sprintf('%s_Series_%d_t%0.3d_c001.tif', prefix, s, t)];
            outputPath = ['../code/experiments/Processed data/' sprintf('%s_Series_%d_t%0.3d_c001.tif', prefix, s, t)];

            cmd = sprintf('convert "%s" -level %.3f%%,%.3f%% "%s"', inputPath, level_black, level_white, outputPath);
            fprintf(fid, [strrep(cmd, '%', '%%') ' 2> /dev/null\n']);

            % Resize the output of BF images
            cmd = sprintf('convert "%s" -resize 25%% "%s"', outputPath, outputPath);
            fprintf(fid, [strrep(cmd, '%', '%%') '\n']);

            cmd = sprintf('convert "%s" -gravity north -stroke none -fill white -pointsize 50 -annotate 0 "%s" "%s"', outputPath, sprintf('BF, T = %.1f h', TT(T(t))+dT(s)), outputPath);
            fprintf(fid, [cmd '\n']);

            cmd = sprintf('convert "%s" -gravity southeast -stroke none -fill white -pointsize 36 -annotate +5+20 "100 µm" "%s"', outputPath, outputPath);
            fprintf(fid, [cmd '\n']);

            cmd = sprintf('composite -compose over -gravity southeast -geometry +5+5 \\\\( -size %dx10 xc:white \\\\) "%s" "%s"', round(100/dp), outputPath, outputPath);
            fprintf(fid, [cmd '\n']);


            % Rescale color and resize output of GFP images
            inputPath = ['Exported data/' sprintf('%s_Series_%d_t%0.3d_c002.tif', prefix, s, t)];
            outputPath = ['../code/experiments/Processed data/' sprintf('%s_Series_%d_t%0.3d_c002.tif', prefix, s, t)];

            level_black = 100*median(BG_GFP(:))/2^16;
            level_white = 100*max(MAX_GFP(:))/2^16;

            cmd = sprintf('convert "%s" -resize 25%% -level %.3f%%,%.3f%% -fill green -tint 100 "%s"', inputPath, level_black, level_white, outputPath);
            fprintf(fid, [strrep(cmd, '%', '%%') ' 2> /dev/null\n']);

            cmd = sprintf('convert "%s" -gravity north -stroke none -fill white -pointsize 50 -annotate 0 "%s" "%s"', outputPath, sprintf('GFP, T = %.1f h', TT(T(t))+dT(s)), outputPath);
            fprintf(fid, [cmd '\n']);

            cmd = sprintf('convert "%s" -gravity southeast -stroke none -fill white -pointsize 36 -annotate +5+20 "100 µm" "%s"', outputPath, outputPath);
            fprintf(fid, [cmd '\n']);

            cmd = sprintf('composite -compose over -gravity southeast -geometry +5+5 \\\\( -size %dx10 xc:white \\\\) "%s" "%s"', round(100/dp), outputPath, outputPath);
            fprintf(fid, [cmd '\n']);

            % Combine into comparrison image
            inputPath =@(c) ['../code/experiments/Processed data/' sprintf('%s_Series_%d_t%0.3d_c00%d.tif', prefix, s, t, c)];
            outputPath = ['../code/experiments/Processed data/' sprintf('B%dX%d_t%0.3d.tif', C(s), sum(C(1:s)==C(s)), T(t))];

            cmd = sprintf('convert +append "%s" "%s" -resize 50%% "%s"', inputPath(1), inputPath(2), outputPath);
            fprintf(fid, [strrep(cmd, '%', '%%') '\n']);

            cmd = sprintf('rm "%s" "%s"', inputPath(1), inputPath(2));
            fprintf(fid, [cmd '\n']);

        end
    end
end

save([dpath '/data.mat'], 'TT', 'dT', 'BG_GFP', 'MAX_GFP', 'radius', 'GFP_radius_0025', 'GFP_radius_0030', 'GFP_radius_0035', 'GFP_radius_0040', 'GFP')

fclose(fid);
