
function salientWatershed_ver2(inputDirectory, fName, outputDirectory)

%****************************
% Parvez Ahammad, JFRC, HHMI
%****************************
% Code to improve watershed via computing salient edges in a hierarchical approach
% step-1: a simple non-local averaging - Buades2005 / CVPR
% step-2: canny with a very low threshold value (low false-negative)
% step-3: enhance image boundaries using distance tf of canny edges
% step-4: watershed partition of the distance map of salient edge image
%
% Notes/Comments
%****************
% 04/13/11: removed the pointwise multiply with original image. This means
% that the watershed will be applied to the saleint edge image at the end.
% 04/14/11: added the changes to remove small edges from the salient edge
% set
% 04/14/11: added interactive setting of the parameters

tic

%%%%%%%%%%
% clear all; clc; close all; pause(0.1)
% inputDirectory = '/Users/ahammadp/Documents/parvez-laptop/work/data/EM_xPrep_segmentation/1000x1000/sub/';
% fName = '01_0.tif';
%inputDirectory = '/Users/ahammadp/Documents/parvez-laptop/work/data/BSDS300/images/train/';
%fName = '41004.jpg';
%%%%%%%%%%%%%%%

%%Processing parameters
%-------------------------------
%NLmeans_flag = input('\n Enter 1 to perform non-local means filtering, 0 otherwise:');
NLmeans_flag=1;
%plot_flag = input('\n Enter 1 to display results, 0 otherwise:');
plot_flag=0;
%salientEdgeThreshold = input('\n Enter the minimum acceptable salient edge length, in pixels:');
salientEdgeThreshold = 3;

texton_filter_support = 25;
NLMeansSigma = 50;
cannyThresh = 1/200;
pbThresh = 1/200;
distMapScale = 2;
%-------------------------------

%% Pre-processing
%add path to access Pb computation routine
addpath('/Users/ahammadp/Documents/parvez-laptop/work/code/segbench/lib/matlab');

if NLmeans_flag == 0
    NLMeansSigma = 0;
end

if ~exist('outputDirectory', 'var')
    outputDirectory = [inputDirectory 'resultsDump/'];
end
if ~exist(outputDirectory, 'dir')
    mkdir([inputDirectory 'resultsDump/']);
end


fprintf('\n Input directory is: %s',inputDirectory);
fprintf('\n Input file name is: %s',fName);
fprintf('\n Output directory is: %s \n\n',outputDirectory);

%im = imcomplement(imread (strcat(inputDirectory,fName)));
im = imread (strcat(inputDirectory,fName));

% Make sure that the input image is 8 bit 2-D image.
if size(im,3)==3
    fprintf('\n Converting RGB image to Grayscale image .. \n');
    im = rgb2gray(imread(strcat(inputDirectory,fName)));
end

%% non-local averaging filter to denoise the image
fprintf('\n Non-local means computation on original image.. \n')
if NLmeans_flag ==1
    filtNLIm = NLmeansfilter(double(im),3,3,NLMeansSigma);
    imwrite(uint8(filtNLIm), [outputDirectory fName '_' 'OriginalNLmeans' '.tif'],'TIFF', 'Compression', 'none');
else
    filtNLIm = im;
end

%% compute canny boundary
fprintf('\n Computing canny edges.. \n')
cannyEdgeIm = bwareaopen(edge(filtNLIm, 'canny', cannyThresh), salientEdgeThreshold);

%% compute probability of boundary
fprintf('\n Computing Pb.. \n')
[pb_orig,~] = pbBGTG2(filtNLIm/max(filtNLIm(:)),'gray');
adjPb = imadjust(histeq(pb_orig));
adjThickPb = pbThicken(adjPb);

%% compute pbEdgeIm using PbThresh
fprintf('\n Combining Pb and Canny.. \n')
pbEdgeIm = adjThickPb>pbThresh;
% trying to implement connected component based stray edge elimination
hybridPbThickCannyEdgeIm = bwareaopen(pbEdgeIm.*cannyEdgeIm, salientEdgeThreshold);


%% enhance original image with estimated boundary map

%old code with pointwise multiply with original image
% enhanceIm = exp(-distHybridIm*distMapScale).*double(imresize(im,size(distHybridIm)));

%new update 04/13/2011
enhanceIm = exp(-bwdist(hybridPbThickCannyEdgeIm)*distMapScale);
enhanceIm2 = exp(-bwdist(cannyEdgeIm)*distMapScale);

fprintf('\n Computing watershed partitions on edge enhanced image.. \n')
ws_enhanceIm = watershed(enhanceIm);
ws_enhanceIm2 = watershed(enhanceIm2);

numSP = max(ws_enhanceIm(:))
numSP2 = max(ws_enhanceIm2(:))

% making composites to overlay edge maps on images
composite1 = plotGreenEdges(im, hybridPbThickCannyEdgeIm, plot_flag);
composite2 = plotRedEdges(im, ws_enhanceIm==0, plot_flag);
composite3 = plotRedGreenEdges(im, ws_enhanceIm==0, hybridPbThickCannyEdgeIm, plot_flag);
composite4 = plotRedGreenEdges(im, watershed(filtNLIm)==0, hybridPbThickCannyEdgeIm, plot_flag);
composite5 = plotRedEdges(im, watershed(filtNLIm)==0, plot_flag);
composite6 = plotRedEdges(im, watershed(filtNLIm)==0, plot_flag);
composite7 = plotRedGreenEdges(im, ws_enhanceIm==0, cannyEdgeIm, plot_flag);
composite8 = plotRedGreenEdges(im, watershed(filtNLIm)==0, cannyEdgeIm, plot_flag);
composite9 = plotRedGreenEdges(im, cannyEdgeIm, hybridPbThickCannyEdgeIm, plot_flag);

composite2_2 = plotRedEdges(im, ws_enhanceIm2==0, plot_flag);
composite3_2 = plotRedGreenEdges(im, ws_enhanceIm2==0, hybridPbThickCannyEdgeIm, plot_flag);
composite7_2 = plotRedGreenEdges(im, ws_enhanceIm2==0, cannyEdgeIm, plot_flag);


% adjThickPb.*(ws_enhanceIm==0) is worth checking!!
if NLmeans_flag == 1
    imwrite(composite1,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'Hybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite2,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsEnhance' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite3,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsEnhanceHybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite4,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsOriginalHybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite5,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsOriginal' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite6,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsOriginalNLFilt' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite7,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsEnhanceCannyEdge' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite8,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsOriginalCannyEdge' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite9,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'CannyEdgeHybrid' '.tif'],'TIFF', 'Compression', 'none');

    imwrite(composite2_2,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsEnhance2' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite3_2,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsEnhance2Hybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite7_2,[outputDirectory fName '_' 'NLMeansSigma' num2str(NLMeansSigma) '_' 'wsEnhance2CannyEdge' '.tif'],'TIFF', 'Compression', 'none');
else
    imwrite(composite1,[outputDirectory fName '_' 'NoNLmeans' '_' 'Hybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite2,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsEnhance' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite3,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsEnhanceHybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite4,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsOriginalHybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite5,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsOriginal' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite6,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsOriginalNLFilt' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite7,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsEnhanceCannyEdge' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite8,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsOriginalCannyEdge' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite9,[outputDirectory fName '_' 'NoNLmeans' '_' 'CannyEdgeHybrid' '.tif'],'TIFF', 'Compression', 'none');

    imwrite(composite2_2,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsEnhance2' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite3_2,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsEnhance2Hybrid' '.tif'],'TIFF', 'Compression', 'none');
    imwrite(composite7_2,[outputDirectory fName '_' 'NoNLmeans' '_' 'wsEnhance2CannyEdge' '.tif'],'TIFF', 'Compression', 'none');
end

%% prepare output structure
%data
result.im = im;
result.filtNLIm = filtNLIm;
result.cannyEdgeIm = cannyEdgeIm;
result.pbEdgeIm = pbEdgeIm;
result.hybridPbThickCannyEdgeIm = hybridPbThickCannyEdgeIm;
result.enhanceIm= enhanceIm;
result.ws_enhanceIm = ws_enhanceIm;
result.enhanceIm2= enhanceIm2;
result.ws_enhanceIm2 = ws_enhanceIm2;

%parameters
result.NLmeans_flag = NLmeans_flag;
result.texton_filter_support = texton_filter_support;
result.NLMeansSigma = NLMeansSigma;
result.cannyThresh = cannyThresh;
result.pbThresh = pbThresh;
result.distMapScale = distMapScale;
result.salientEdgeThreshold = salientEdgeThreshold;

eval(sprintf('save %s result',strcat(outputDirectory,strcat(fName,'_Result.mat'))));
%clear composite1 composite2 composite3 composite4 enhanceIm distHybridIm
%return
toc