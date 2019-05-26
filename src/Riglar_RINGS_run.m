function Riglar_RINGS_run(varargin)

%% parse inputs

ip = inputParser;
ip.addParameter('fitSingleColonyMethod', @Riglar_RINGS_fit);
ip.KeepUnmatched = true;

ip.parse(varargin{:});
p = ip.Results;
pp = ip.Unmatched;

%%

% get list of folders with images, corresponding to plates
inputFolderName = uigetdir('','Select parent folder containing YFP and CFP folders');

yfp_path = [inputFolderName filesep 'YFP'];
cfp_path = [inputFolderName filesep 'CFP'];

platesList = dir(yfp_path); % use yfp as primary path.  assumes cfp has the same plates (with the same names)

% remove non-directories
dirCheck1 = vertcat(platesList.isdir);
dirCheck2 = not(strcmp({platesList.name}, '.') | strcmp({platesList.name}, '..'));
dirCheck  = logical(dirCheck1(:) .* dirCheck2(:));
platesList = platesList(dirCheck);

for plateID = 1:numel(platesList)
    
    fnameList = dir([yfp_path filesep platesList(plateID).name]);

    % remove directories
    imFileCheck1 = not(vertcat(fnameList.isdir));
    
    % remove non-image files
    imFileCheck2 = zeros(size(fnameList));
    for fileID = 1:numel(fnameList)
        imFileCheck2(fileID) = any(~cellfun(@isempty, regexp(fnameList(fileID).name, {'.tif','.tiff','.jpg','.jpeg','.png'}))); % simply add other file extensions here
    end
    
    imFileCheck = logical(imFileCheck1(:) .* imFileCheck2(:));
    fnameList = fnameList(imFileCheck);
    
    yfp_fnames = {fnameList.name};
    for i = 1:numel(yfp_fnames)
        cfp_fnames{i} = strrep(yfp_fnames{i},'YFP','CFP');
    end
    
    for fileID = 1:numel(fnameList)
                
        yfp_im = double(imread([yfp_path filesep platesList(plateID).name filesep yfp_fnames{fileID}]));
        cfp_im = double(imread([cfp_path filesep platesList(plateID).name filesep cfp_fnames{fileID}]));
        im(:,:,1) = yfp_im;
        im(:,:,2) = cfp_im;
        
        fitInfo{fileID,plateID} = p.fitSingleColonyMethod(im, pp);
        
    end
end

%% write to table for output

% initialize table for output
output = array2table(NaN(size(fitInfo)));
for i = 1:numel(platesList)
    output.Properties.VariableNames{i} = ['t_' platesList(i).name];
end

for fileID = 1:size(fitInfo,1)
    for plateID = 1:size(fitInfo,2)
        
        if not(isempty(fitInfo{fileID,plateID})) && (fitInfo{fileID,plateID}.qualityFlag == 0)
            output{fileID,plateID} = mod(fitInfo{fileID,plateID}.theta_0, 2*pi);
        else
            output{fileID,plateID} = NaN;
        end
    end
end

%% visualize result
    
cf = figure;
hold on,

for plateID = 1:numel(platesList)
    
    scatter(ones(1,numel(output.(plateID)))*(plateID-1), output.(plateID))
    
end

axis([-0.5 (numel(platesList)-0.5) 0 2*pi])
cf.CurrentAxes.XTick = 0:numel(platesList)-1;
cf.CurrentAxes.XTickLabel = {platesList.name};
cf.CurrentAxes.XTickLabelRotation = 45;

%% save output

outputFolderName = uigetdir('','Select folder to write calibration data');

writetable(output, [outputFolderName filesep 'theta_0.csv'])