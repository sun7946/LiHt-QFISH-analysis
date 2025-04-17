classdef dot_TimelessM
    %%
    %  Description:
    %             imgage processing class for static images, see readme_TimelessM.md for details
    %             also for data visual and trouble shooting             
    %
    %  Modified:
    %             created by: Min Lab, Oct 11th, 2022
    %
    %  Ver Tag:
    %             V1
    %
    %  Copyright:
    %             copyright statment goes here
    %

    %% Properties
    % Public, non-tunable properties
    properties(SetAccess = public)
        path;
        version = '1.0';
        updatedTime = '2022-10-11';
    end

    properties (SetAccess = public, GetAccess = public)
        imagePath;
        maskPath;
        outPath;
        fileName;

        row;
        col;
        site;
        shot;
        fileType = 'nd2';
        channel;

        % image parameters
        height;
        width;
        bias = 1;
        cmos = 100;

        % data file
        dataFile;
%         cellID;
    end

    %% Methods
    methods (Access = public)
        %% single image preprocessing: reading, shading correction and background removal
        function obj = initiate(obj, varargin)
            tic;
            disp(['The TimelessM version is ', obj.version, ', updated on ', obj.updatedTime, '.']);

            %% initialize parser
            p = inputParser;

            % add the function name to the input parser
            p.FunctionName = mfilename;

            %% default tag/value pairs
            p.addParameter('imagePath',         obj.imagePath,      @ischar);
            p.addParameter('maskPath',          obj.maskPath,       @ischar);
            p.addParameter('outPath',           obj.outPath,        @ischar);
            p.addParameter('fileType',          obj.fileType,       @ischar);
            p.addParameter('fileName',          obj.fileName,       @ischar);
            p.addParameter('dataFile',          obj.dataFile,       @ischar);

            p.addParameter('row',               obj.row,            @isnumeric);
            p.addParameter('col',               obj.col,            @isnumeric);
            p.addParameter('site',              obj.site,           @isnumeric);
            p.addParameter('channel',           obj.channel,        @isnumeric);
%             p.addParameter('cellID',            obj.cellID,         @isnumeric);

            p.addParameter('bias',              obj.bias,           @isnumeric);
            p.addParameter('cmos',              obj.cmos,           @isnumeric);

            % parse the variable argument list
            p.parse(varargin{:});

            obj.imagePath       = p.Results.imagePath;
            obj.maskPath        = p.Results.maskPath;
            obj.outPath         = p.Results.outPath;
            obj.bias            = p.Results.bias;
            obj.cmos            = p.Results.cmos;
            obj.channel         = p.Results.channel;
            obj.site            = p.Results.site;
            obj.fileType        = p.Results.fileType;
            obj.fileName        = p.Results.fileName;
            obj.dataFile        = p.Results.dataFile;
            obj.row             = p.Results.row;
            obj.col             = p.Results.col;

            %% define output folder
            if isempty(obj.outPath)
                obj.outPath = pwd;
            end

            if ~exist(obj.outPath,'dir')
                mkdir(obj.outPath);
            end

            %% parse data file
            if ~isempty(obj.dataFile)
                [~,fn] = fileparts(obj.dataFile);
                fileParts = split(fn, '_');
                obj.row = str2double(fileParts{1});
                obj.col = str2double(fileParts{2});
                obj.site = str2double(fileParts{3});
            end

            %% parse row and col number
            if isempty(obj.row) || isempty(obj.col)
                if isempty(obj.fileName)
                    error('No info of input file!');
                else
                    if strcmp(obj.fileType, 'nd2') 
                        %Find a capital letter, followed by two digits
                        startIdx = regexp(obj.fileName,'[A-Z][0-9][0-9]');
                    
                        %Convert the first letter into a number
                        obj.row = int8(obj.fileName(startIdx)) - 64; %'A' = 65
                    
                        %Convert the next two digits into the column number
                        obj.col = str2double(fn(startIdx+1:startIdx+2));
                    else
                        newStr = split(obj.fileName,'_');
                        obj.row = str2double(newStr(1));
                        obj.col = str2double(newStr(2));
                    end
                end
            end
            
            %% set shot
            obj.shot = [num2str(obj.row),'_',num2str(obj.col),'_',num2str(obj.site)];

            %% get file name
            if isempty(obj.fileName) && strcmp(obj.fileType, 'nd2')
                myAlphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
                rowLetter = myAlphabet(obj.row);
                obj.fileName = sprintf([rowLetter,'%02d'], obj.col);

                currFiles = dir([fullfile(obj.imagePath, filesep),'*nd2']);

                ndFileNames = {currFiles.name};
                find_match = regexp(ndFileNames, obj.fileName,'match');
                find_match = ~cellfun('isempty',find_match);
                file_match = find(find_match>0,1);

                if ~isempty(file_match)
                    obj.fileName = currFiles(file_match).name;
                else
                    error('No files found! Please check image path and file name.')
                end
            end

            %% get the height and width of the images
            dummyImage = readImage(obj,1,obj.channel); 
            [obj.height, obj.width] = size(dummyImage);    

            %% get site numbers
            if isempty(obj.site)
                if strcmp(obj.fileType, 'nd2')
                    myBfReader = BioformatsImage(fullfile(obj.imagePath, obj.fileName));
                    obj.site = 1:myBfReader.seriesCount;
                else
                    tifFile = dir([fullfile(obj.maskPath, filesep),'*tif']);
                    testFlag = 1;
                    testSite = 1;
                    testFile = 1;
                    while testFlag
                        tifName = tifFile(testFile).name;
                        tifNameParts = split(tifName,'_');
                        if testSite<tifNameParts(3)
                            testSite = tifNameParts(3);
                            testFile = testFile + 1;
                        else
                            obj.site = 1:testSite;
                            testFlag = 0;
                        end
                    end
                end
            end
            
            %% get channel numbers
            if isempty(obj.channel)
                if strcmp(obj.fileType, 'nd2')
                    myBfReader = BioformatsImage(fullfile(obj.imagePath, obj.fileName));
                    obj.channel = 1:myBfReader.sizeC;
                else
                    tifFile = dir([fullfile(obj.imagePath, filesep),'*tif']);
                    testFlag = 1;
                    testChannel = 1;
                    testFile = 1;
                    while testFlag
                        tifName = tifFile(testFile).name;
                        tifNameParts = split(tifName,'_');
                        if testSite<tifNameParts(4)
                            testChannel = tifNameParts(4);
                            testFile = testFile + 1;
                        else
                            obj.channge = 1:testChannel;
                            testFlag = 0;
                        end
                    end
                end
            end
            toc;
        end

        %% preprocessing - export backgound corrected tif for segmentation etc.
        function obj = prePro(obj, bgRemoveMethod)
            if nargin > 1 && ismember(bgRemoveMethod, ["global", "local", "rlocal"])
                bgMethod = bgRemoveMethod;
            else
                bgMethod = 'global';
            end

            for i = obj.site
                fprintf('frame %0.0f\n',i);
                rawImg = readImage(obj, i, obj.channel(1));
                correctedImg = (rawImg - obj.cmos)./obj.bias(:,:,obj.channel(1));
                realImg = removeBackground(obj, correctedImg, bgMethod);
                imwrite(uint16(realImg), fullfile(obj.outPath, [obj.shot,'_',num2str(obj.channel(1)),'_',num2str(i),'.tif']));
            end

        end

        %% extract signal
        function obj = IFSignalExtract(obj)
            %% identify the naming pattern of mask files
            maskFile = dir([fullfile(obj.maskPath, filesep),'*tif']);
            dummyName = maskFile(1).name;
            maskFileParts = split(dummyName,'_');
            startFlag = 1;
            if str2double(maskFileParts{3}) == 0
                startFlag = 0;
            end
            sepIndex = regexp(dummyName,'_');
            maskEnd = dummyName(sepIndex(3):end);

            %% process each file of view
            se = strel('disk',50); 
            real = NaN(obj.height, obj.width, length(obj.channel));

            for i = obj.site

                shot = [num2str(obj.row),'_',num2str(obj.col),'_',num2str(i)];
                
                % read and correct image, remove background
                for j = 1:length(obj.channel)
                    raw = readImage(obj, i, obj.channel(j));
                    corrected = (raw - obj.cmos)./obj.bias(:,:,j);
                    real(:,:,j) = imtophat(corrected,se);
                end

                % load mask
                nuc_mask = single(imread([obj.maskPath, [num2str(obj.row),'_',num2str(obj.col),'_',num2str(i-1+startFlag)] ,maskEnd])); 
                
                % basic nuclear info
                nuc_info = regionprops(nuc_mask,'Area','Centroid','PixelIdxList');
                ringWidth = 3;
                ringLabel = getCytoring(obj, nuc_mask, ringWidth);
                ring_info = regionprops(ringLabel,'PixelIdxList');
                
                % extract cell info (this part take most of the processing time)
                numcells = length(nuc_info);   

                parameterNum = length(obj.channel)*4+3;
                cellinfo = NaN(numcells,  parameterNum);
             
                for cc = 1:numcells
                    cellinfo(cc,1:2)    = nuc_info(cc).Centroid;
                    cellinfo(cc,3)         = nuc_info(cc).Area;

                    for k = 1:length(obj.channel)
                        currChannel = obj.channel(k);
                        currReal = real(:,:,currChannel);
                        cellinfo(cc, 4*k)           = median(currReal(nuc_info(cc).PixelIdxList));
                        cellinfo(cc, 4*k+1)         = median(currReal(ring_info(cc).PixelIdxList));
                        cellinfo(cc, 4*k+2)         = sum(currReal(nuc_info(cc).PixelIdxList));
                        cellinfo(cc, 4*k+3)         = mean(currReal(nuc_info(cc).PixelIdxList));
                    end
                end

                % convert matrix to struct
                a = [compose('ch%d_median', 1:length(obj.channel)); compose('ch%d_ring_median', 1:length(obj.channel)); compose('ch%d_sum', 1:length(obj.channel)); compose('ch%d_mean', 1:length(obj.channel))];
                
                structName = ['x';'y';'area'; a(:)];
                cells = cell2struct(num2cell(cellinfo, parameterNum), structName, 2);
        
                % store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isempty(nuc_info)
                    save([obj.outPath,shot,'_IF.mat'],'cells');
                end
            end
        end

        %% match cells of interest onto images.
        function obj = imgDataMatch(obj, idList)
            h = load(obj.dataFile);
            blank = zeros(obj.height, obj.width);

            if ~isempty(obj.maskPath)

                % identify the naming pattern of mask files
                maskFile = dir([fullfile(obj.maskPath, filesep),'*tif']);
                dummyName = maskFile(1).name;
                maskFileParts = split(dummyName,'_');
                startFlag = 1;
                if str2double(maskFileParts{3}) == 0
                    startFlag = 0;
                end
                sepIndex = regexp(dummyName,'_');
                maskEnd = dummyName(sepIndex(3):end);

                 % load mask
                nuc_mask = single(imread([obj.maskPath, [num2str(obj.row),'_',num2str(obj.col),'_',num2str(obj.site-1+startFlag)] ,maskEnd])); 
                nuc_info = regionprops(nuc_mask,'Area','Centroid','PixelIdxList');
                for i = 1:length(idList)
                    cellID = idList(i);
                    blank(nuc_info(cellID).PixelIdxList) = cellID;
                end
                figure, imshow(blank,[]);
                imwrite(uint16(blank),'test1.tif');
            else
                try
                    coor = cat(1, h.cells.xy_coordinates);
                    x = floor(coor(:,2));
                    y = floor(coor(:,1));
                catch
                    x = floor(cat(1, h.cells.y));
                    y = floor(cat(1, h.cells.x));
                end
                
                for i = 1:length(idList)
                    cellID = idList(i);
                    blank(x(cellID),y(cellID)) = 1;
                end
                test = imfilter(blank,fspecial('disk',6),'symmetric');
                figure, imshow(test,[]);
                imwrite(uint16(test*1000),'test1.tif');
            end

            
        end

        %% remove background
        function [real, bg] = removeBackground(~, inputImg,bgRemoveMethod)
            % estimate the average background
            [density, intensity] = ksdensity(inputImg(:));
            [~, peakIndex] = findpeaks(density,'SortStr','descend'); 
            bg = intensity(peakIndex(1));

            switch bgRemoveMethod

                case 'global'
                    background = prctile(inputImg(:), 20); % default parameter: background is set as the bottom 10 percentile
                    real = inputImg - background;

                case 'rglobal'
                    % subtract the average background from all pixels in the image
                    real = inputImg - bg;

                case 'local'
                    se      = strel('disk', 50); % default parameter: diameter of image opening is set to 50
                    real    = imtophat(inputImg,se);

                case 'rlocal'
                    % do imtophat (aka 'local') first
                    se          = strel('disk', 50); % default parameter: diameter of image opening is set to 50
                    realtemp    = imtophat(inputImg,se);

                    % subtract the peak intensity from all pixels in the image
                    [density, intensity] = ksdensity(realtemp(:));
                    [~, peakIndex] = findpeaks(density,'SortStr','descend');
                    toSubtract = intensity(peakIndex(1));
                    if toSubtract<2 % in case some pixel intensity stag at 0 (extremely rare after imtophat), use the second peak
                        toSubtract = intensity(peakIndex(2));
                    end
                    real = realtemp - toSubtract;
            end

        end


        %% read images
        function raw = readImage(obj,site,channelNum)
            if strcmp(obj.fileType,'nd2')

                currentFile = fullfile(obj.imagePath, obj.fileName);
                myBfReader = BioformatsImage(currentFile);

                raw = double(myBfReader.getXYplane(channelNum,site, 1));

            else
                try
                    raw = double(imread(fullfile(obj.imagePath,sprintf('%s_%s_%s.tif', obj.shot, num2str(channelNum), num2str(site)))));
                catch
                    raw = double(imread(fullfile(obj.imagePath,sprintf('%s_%s_%s.tif', obj.shot, num2str(site), num2str(channelNum)))));
                end
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% calculate MIP
        function  obj = MIP(obj)
            currentFile = fullfile(obj.imagePath, obj.fileName);
            reader = BioformatsImage(currentFile);
            % Obtain dimensional information of images
            numZ = reader.sizeZ();
            numC = reader.sizeC();
            numT = reader.sizeT();
            MIP = zeros(reader.height(), reader.width(), 'double');  % Create a matrix for storing the maximum intensity projection
            % Loop through the Z-stack and calculate the maximum intensity projection
            for z = 1:numZ
                plane= double(reader.getPlane([z, obj.channel, numT, obj.site]));  % Reading image data
                MIP = max(MIP, plane);  % MIP
            end
            imwrite(uint16(MIP), fullfile(obj.outPath, [obj.shot,'_',num2str(obj.channel),'_',num2str(numT),'.tif']));  % save MIP
        end

        %% autofocus
        function bestImage = Autofocus(obj,channel)
            currentFile = fullfile(obj.imagePath, obj.fileName);
            reader = BioformatsImage(currentFile);
            % Obtain dimensional information of images
            numZ = reader.sizeZ();
            numC = reader.sizeC();
            numT = reader.sizeT();
            % initialization
            bestSlice = 1;
            bestFocusScore = 0;
            % Traverse Z-axis slices and calculate focus score
            for slice = 1:numZ
                image = double(reader.getPlane([slice, channel, numT, obj.site]));  % read image
                %  Calculate focus score: Sobel gradient
                gradientImage = imgradient(image);
                focusScore = sum(gradientImage(:));
                % Update the clearest image information
                if focusScore > bestFocusScore
                    bestFocusScore = focusScore;
                    bestSlice = slice;
                end
            end
            bestImage = double(reader.getPlane([bestSlice, channel, numT, obj.site]));  % Obtain the clearest image slice
%             imwrite(uint16(bestImage), fullfile(obj.outPath, [obj.shot,'_',num2str(obj.channel),'_',num2str(numT),'.tif']));  % save image
        end
        %%

        %% dot signal processing
        function obj = dot_process(obj,bg_parameter,threshold_seg,mask_outputFolder)
            img_name = [obj.shot,'_',num2str(obj.channel),'_1.tif'];
            raw_img = single(imread([obj.imagePath,img_name]));
%             raw = raw_img;
%             raw(raw>2000) = 2000;
            %    figure; imshow(raw,[]);impixelinfo
            %% Remove the internal background of cells
            se_open = strel('disk', bg_parameter);   % Remove the parameters of intracellular background, the smaller the value, the more stringent the removal of background
            Rst_open = imopen(raw_img, se_open);    % imopen: corrode first and then expand
            %    figure; imshow(Rst_open,[]);impixelinfo
            img1 = uint16(raw_img - Rst_open);    %% folder：dilation-FITC_imopen
            img1(Rst_open<0) = 0;  % Assign pixel values less than 0 to 0
            %    figure; imshow(img1,[]);impixelinfo
            imwrite(uint16(img1), fullfile(obj.outPath, img_name));
            %% filter image - alternative method: imtophat
            se = strel('disk', 5); 
            filtered_img = imtophat(img1, se);
            %     imwrite(uint16(filtered_img), fullfile(outputFolder, [img_name,'_2_202','.tif']));
            %% mask the dots
            dots_mask = imbinarize(filtered_img,threshold_seg);  % Segmentation threshold, the smaller the value, the lower the threshold
            mask = bwlabel(dots_mask,4);
            imwrite(uint16(mask),fullfile(mask_outputFolder, img_name));
            toc
        end
        
        %% extract signal: telomere
        function obj = telomere_SignalExtract(obj,telomerePath,telomereMaskPath,telomere_channel)
            %% identify the naming pattern of mask files
            maskFile = dir([fullfile(obj.maskPath, filesep),'*191.tif']);
            dummyName = maskFile(1).name;
            maskFileParts = split(dummyName,'_');
            startFlag = 1;
            if str2double(maskFileParts{3}) == 0
                startFlag = 0;
            end
            sepIndex = regexp(dummyName,'_');
            maskEnd = dummyName(sepIndex(3):end);

            %% process each file of view
            se = strel('disk',50); 
            real = NaN(obj.height, obj.width, length(obj.channel));

            for i = obj.site

                shot = [num2str(obj.row),'_',num2str(obj.col),'_',num2str(i)];
                channel_length = 1:length(obj.channel);
                isTelomere = ismember(channel_length,telomere_channel);
                channel_ex_tel = channel_length(~isTelomere);

                % read and correct image, remove background
                for j = channel_ex_tel(1):channel_ex_tel(length(channel_ex_tel))
                    raw = Autofocus(obj,j);
                    corrected = (raw - obj.cmos)./obj.bias(:,:,j);
                    real(:,:,j) = imtophat(corrected,se);
                end

                % load mask
                nuc_mask = single(imread([obj.maskPath, [num2str(obj.row),'_',num2str(obj.col),'_',num2str(i-1+startFlag)] ,maskEnd])); 
                
                % load telomere
                telomere = single(imread([telomerePath, obj.shot,'_',num2str(telomere_channel),'_191.tif']));
                real(:,:,telomere_channel) = telomere;

                % load telomere mask
                telomere_mask = single(imread([telomereMaskPath, obj.shot,'_',num2str(telomere_channel),'_191.tif']));

                % basic nuclear info
                nuc_info = regionprops(nuc_mask,'Area','Centroid','PixelIdxList','PixelList','ConvexHull');
                ringWidth = 3;
                ringLabel = getCytoring(obj, nuc_mask, ringWidth);
                ring_info = regionprops(ringLabel,'PixelIdxList');
                ch_telomere_value = regionprops(telomere_mask,telomere,'PixelValues','PixelIdxList','Centroid','Area','PixelList');
                
                % extract cell info (this part take most of the processing time)
                numcells = length(nuc_info);   

                parameterNum = (length(obj.channel)-1)*4+13;
                cellinfo = NaN(numcells,  parameterNum);
             
                for cc = 1:numcells
                    cellinfo(cc,1:2)    = nuc_info(cc).Centroid;
                    cellinfo(cc,3)      = nuc_info(cc).Area;

                    for k = 1:length(obj.channel)-1
                        currChannel = obj.channel(k);
                        currReal = real(:,:,currChannel);
                        cellinfo(cc, 4*k)      = median(currReal(nuc_info(cc).PixelIdxList));
                        cellinfo(cc, 4*k+1)    = median(currReal(ring_info(cc).PixelIdxList));
                        cellinfo(cc, 4*k+2)    = sum(currReal(nuc_info(cc).PixelIdxList));
                        cellinfo(cc, 4*k+3)    = mean(currReal(nuc_info(cc).PixelIdxList));
                    end
                end

                % convert matrix to struct
                a = [compose('ch%d_median', channel_ex_tel(1):channel_ex_tel(length(channel_ex_tel))); ...
                     compose('ch%d_ring_median', channel_ex_tel(1):channel_ex_tel(length(channel_ex_tel))); ...
                     compose('ch%d_sum', channel_ex_tel(1):channel_ex_tel(length(channel_ex_tel))); ...
                     compose('ch%d_mean', channel_ex_tel(1):channel_ex_tel(length(channel_ex_tel)))];

                structName = ['x';'y';'area'; a(:);...
                              'ch_telomere_mean';...
                              'ch_telomere_median';...
                              'ch_telomere_sum';...
                              'ch_telomere_dot_Centroid';...
                              'ch_telomere_dot_area';...
                              'ch_telomere_dot_median_in_median';...
                              'ch_telomere_dot_median_in_sum';...
                              'ch_telomere_dot_sum_in_median';...
                              'ch_telomere_dot_sum_in_sum';...
                              'ch_telomere_dot_number';];
                cells = cell2struct(num2cell(cellinfo, parameterNum), structName, 2);
                tic
                for cc = 1:numcells
                    
                    num = 1;
                    dot_xy = [];
                    dot_xy = cat(1,ch_telomere_value.Centroid);
                    [in, on] = inpolygon(dot_xy(:,1), dot_xy(:,2), nuc_info(cc).ConvexHull(:,1), nuc_info(cc).ConvexHull(:,2));
                    loc = in | on;
                    dot_index = find(loc==1);
                    %Determine whether FISH point Centroid intersects with the cell nucleus mask,
                    % that is, whether the number of rows and columns in the intersection result is 0. 
                    % If there is an intersection, return 1
                    for j = 1 : length(dot_index)
                            cells(cc).ch_telomere_mean(num,1)           = mean(ch_telomere_value(dot_index(j)).PixelValues);
                            cells(cc).ch_telomere_median(num,1)         = median(ch_telomere_value(dot_index(j)).PixelValues);
                            cells(cc).ch_telomere_sum(num,1)            = sum(ch_telomere_value(dot_index(j)).PixelValues);
                            cells(cc).ch_telomere_dot_Centroid(num,1:2) = ch_telomere_value(dot_index(j)).Centroid;
                            cells(cc).ch_telomere_dot_area(num,1)       = ch_telomere_value(dot_index(j)).Area;
                            num = num+1;
                    end
                    cells(cc).ch_telomere_dot_median_in_median = median(double(cells(cc).ch_telomere_median));
                    cells(cc).ch_telomere_dot_median_in_sum    = median(double(cells(cc).ch_telomere_sum));
                    cells(cc).ch_telomere_dot_sum_in_median    = sum(double(cells(cc).ch_telomere_median));
                    cells(cc).ch_telomere_dot_sum_in_sum       = sum(double(cells(cc).ch_telomere_sum));
                    cells(cc).ch_telomere_dot_number           = length(cells(cc).ch_telomere_median);
                end
                toc
                % store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if ~isempty(nuc_info)
                    save([obj.outPath,obj.shot,'_telomere.mat'],'cells');
                end
            end                            
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

    end

    methods (Access = private)

        %% get cytoring label from nuc_label
        function ring_label = getCytoring(~, nuc_label,ringwidth)
            %%% define cytosolic ring radii %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            nuc_label = single(nuc_label);
            nuc_mask=nuc_label>0;

            %%% define midband of cytoring with same labels as nuclei %%%%%%%%%%%%%%%%%
            ringedgeouter = imdilate(nuc_label,strel('disk',ringwidth,0));  %outerrad one pixel greater than midrad
            outer_mask = ringedgeouter>0;
            ring_mask = outer_mask - nuc_mask;
            cytoring = ringedgeouter.*ring_mask;       %define cytoring with label ID  %exclude those pixels that are part of nuclear of one cell and cyto_ring of another, MM20170410

            %%% detect cell-cell borders and inflate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            borders = bwmorph(nuc_label,'bothat');
            borders = imdilate(borders,strel('disk',ringwidth,0)); % change 2 to ringwidth - exclude all overlaping pixels, MM20170410

            %%% clarify boundaries %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ring_label=cytoring.*~borders;

            %%%%%%%% reconstitute absent rings %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cru=unique(nuc_label);
            fcru=unique(ring_label);
            noring=setdiff(cru,fcru);

            if ~isempty(noring)
                for i=noring'
                    ringpixel = logical((ringedgeouter==i)-(nuc_label==i));
                    ring_label(ringpixel)=i;
                end
            end
        end

    end

end
