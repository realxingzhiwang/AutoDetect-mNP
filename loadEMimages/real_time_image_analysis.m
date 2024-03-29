% Please do not modify these codes unless you are doing TEM experiments
figure
ButtonHandle = uicontrol('Style', 'PushButton', ...
                         'String', 'Stop loop', ...
                         'Callback', 'delete(gcbf)');

file = 'Y:\Wang\RealTimeAnalysis';
dir_content = dir(file);
filenames = {dir_content.name};
current_files = filenames;
disp(['Monitoring ' file]);

while true
    dir_content = dir(file);
    filenames = {dir_content.name};
    new_files = setdiff(filenames,current_files);
    if ~isempty(new_files)
        disp(['Reading file ' new_files{:}]);
        current_files = filenames;
        pause(1);
        path = fullfile(file, new_files{:});
        [image, image_bw, features] = loadEMimages(path);
        figure('NumberTitle', 'off', 'Name', new_files{:})
        subplot(1,2,1)
        imshow(image)
        title('Original image')
        subplot(1,2,2)
        imshow(image_bw)
        title('Binary image')
        disp('Reading complete')
    end

    if ~ishandle(ButtonHandle)
        disp('Monitoring ended');
        break;
    end
    pause(0.01);


end