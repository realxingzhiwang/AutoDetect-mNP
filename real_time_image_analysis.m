% Please do not modify these codes unless you are doing TEM experiments
addpath(genpath('Hu Moments'), 'iterativeclustering',genpath('loadEMimages'),'naivebayes','ruecs')

close all
figure
ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'delete(gcbf)');

file = '\\tsclient\E\User Data\Xingzhi Wang\20241005\mixture-grid-search2\BW';
dir_content = dir(file);
filenames = {dir_content.name};
current_files = filenames;
disp(['Monitoring ' file]);
% features_org = [];
% particles = {};
% moments = [];
colors = [68 133 255;
        133 255 68;
        255 68 133;
        255 157 37;
        157 37 255;
        37 255 157;
        176 176 176
        117 138 155]/255*0.9;
figure

while true
    dir_content = dir(file);
    filenames = {dir_content.name};
    new_files = setdiff(filenames,current_files);
    if ~isempty(new_files)
        for i = length(new_files)
            if all(new_files{i}(end-2:end)=='tif')
                disp(['Reading file ' new_files{i}]);
                current_files = filenames;
                pause(1);
                path = fullfile(file, new_files{i});
                image_loading = @loadtiff; %@ReadDMFile for dm4, @loadtiff for other formats
                image_segmentation = @identity; %@imagekmeans for mNPs, @combinedthresh for QDs, @identity if inputs are binary images
                area_threshold = 100;
                [image, image_bw, features, particles_ite, unit, ~, moments_ite] = loadEMimages(path, image_loading, image_segmentation,area_threshold);
                features_org = [features_org; features];
                particles = [particles particles_ite{1}];
                moments = [moments; moments_ite];
                % figure('NumberTitle', 'off', 'Name', new_files{i})
                % subplot(1,2,1)
                % imshow(image)
                % title('Original image')
                % subplot(1,2,2)
                % imshow(image_bw)
                % title('Binary image')
                disp('Reading complete')
            end
        end

        if length(particles)>100
            summary_f = core_algorithm(features_org,particles,moments,area_threshold, unit);
            classes_f = summary_f.classification;
            features_f = summary_f.features;
            particles_f = summary_f.particle_shapes;
        
        
            %colored_particles_all = {};
            for j = 1:max(classes_f)
                particles_plot_f = particles_f(classes_f==j);
                colored_particles_f = cell(size(particles_plot_f));
                for i = 1:length(colored_particles_f)
                    colored_particles_f{i} = cat(3, particles_plot_f{i}*colors(j, 1),...
                        particles_plot_f{i}*colors(j, 2), particles_plot_f{i}*colors(j, 3))*1.5;
                end
    
                subplot(1,max(classes_f),j)

                if length(colored_particles_f) < 100 
                    montage(colored_particles_f(1:end), 'BorderSize', [1 1], 'ThumbnailSize', [1 1]*256)
                else
                    subset = mod(1:length(colored_particles_f), 10)==0;
                    montage(colored_particles_f(subset), 'BorderSize', [1 1], 'ThumbnailSize', [1 1]*256)
                    title('x10')
                end
                %colored_particles_all = [colored_particles_all colored_particles_f];
            end
        end
        
    end

    if ~ishandle(ButtonHandle)
        disp('Monitoring ended');
        break;
    end
    pause(0.01);


end