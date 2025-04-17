clearvars
%% CMOS of the camera
h = load('\\172.19.245.39\minlab_data\Groupfolder\Microscope\Gryffindor\CMOS\cmos.mat');
cmos = h.cmos;
width = size(cmos,1); 
heigth = size(cmos,2);
biasAll = ones(width,heigth,4);

%% bias of interest channel
h = load('\\172.19.245.39\minlab_data\Groupfolder\Microscope\Gryffindor\illumination_bias\Gryf_DAPI-20x-Fluorescein_2uM.mat');
biasAll(:,:,1) = h.bias;
h = load('\\172.19.245.39\minlab_data\Groupfolder\Microscope\Gryffindor\illumination_bias\Gryf_FITC_20x_211122.mat');
biasAll(:,:,2) = h.bias;
% h = load('\\192.168.105.249\minlab_data\Groupfolder\Microscope\Gryffindor\illumination_bias\Gryf_FITC_20x_211122.mat');
% biasAll(:,:,3) = h.bias;

%% creat combinations (equivalent to combvec function in the Deep Learning Toolbox)
nd2Path = '\\172.19.245.39\minlab_image1\Qingyang\20250225_telomere_timelapse\fixed_telomere\20250228_131311_529\';
MIP_outPath = 'E:\test\';
bg_prepro_outPath = 'E:\test\';
imopen_outpPath = 'E:\test\';
telomere_mask_outPath = 'E:\test\';
row = 8;
col = 10;
site = 1;
channel = 2;    % telomere channel
bg_parameter = 30;       % 20x for basal, Cellular background parameters
threshold_seg = 0.009;    % 20x for basal，Segmentation threshold
% bg_parameter = 20;       % 20x for basal, Cellular background parameters
% threshold_seg = 0.01;    % 20x for basal，Segmentation threshold
d = combvec_mm(row,col,site);
%% MIP-maximum intensity projection
for i = 1:size(d,1)
    row = d(i,1);
    col = d(i,2);
    site = d(i,3);
    shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
    obj = dot_TimelessM;
    obj = initiate(obj,'row', row,...
        'col',col,...
        'site',site,...
        'channel',channel,...
        'imagePath',nd2Path,...
        'outPath', MIP_outPath);
    obj = MIP(obj);
end
%% Remove medium background
% p = parpool(6);
for i = 1:size(d,1)
    row = d(i,1);
    col = d(i,2);
    site = d(i,3);
    obj = dot_TimelessM;
    obj = initiate(obj,'row', row,...
        'col',col,...
        'site',site,...
        'channel',channel,...
        'bias',biasAll,...
        'cmos', cmos,...
        'fileType','tif',...
        'imagePath',MIP_outPath,...
        'outPath', bg_prepro_outPath);
    obj = prePro(obj,"rlocal");
end
% delete(p)
%% Remove cellular background & achieve telomere mask
for i = 1:size(d,1)
    row = d(i,1);
    col = d(i,2);
    site = d(i,3);
    obj = dot_TimelessM;
    obj = initiate(obj,'row', row,...
        'col',col,...
        'site',site,...
        'channel',channel,...
        'fileType','tif',...
        'imagePath',bg_prepro_outPath,...
        'outPath', imopen_outpPath);
    obj = dot_process(obj,bg_parameter,threshold_seg,telomere_mask_outPath);
end