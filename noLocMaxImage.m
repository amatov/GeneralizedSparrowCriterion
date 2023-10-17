function noLocMaxImage(nbImages,data_which)

% noLocMaxImage substracts from .TIF files GKs 
% positions of the primary and secondary speckles
% and produces the resulting .TIF images
%
% SYNOPSIS   noLocMaxImage(nbImages,data_which)
%            if there is no input parameter - debug mode
%
% INPUT      nbImages   : a number of images to be processed
%            data_which : noise parameters selection
%                         1 - spindle meta red wt 10s 8bit
%                         2 - 158a actin
%                         3 - actin retro
%
% OUTPUT     none       :  
% 
% REMARKS    saves the resulting iamges on disc into
%            a user selected directory
%
% DEPENDENCES   noLocMaxImage uses { }
%               noLocMaxImage is used by { }
%
% Alexandre Matov, November 20th, 2002

DEBUG=0;

if nargin==0
    DEBUG=1;
    nbImages=1;
    data_which=4;
end

s=length(num2str(nbImages));
strg=sprintf('%%.%dd',s); 

SIG=1.88;
enhTriang=0;
userROIbw = [];

switch data_which
    case 1 
        noiseParam=[1.96/1.22 0.02000 2e-4 0.0981 1.96]; % Meta Red Spindle WT 10s 8bit
    case 2
        noiseParam=[1.96/3.55 5.4103e-004 2e-4 0.0295 1.96]; % 158a
    case 3
        noiseParam=[1.96/2.95230 0.03656 1e-4 0.14540 1.96]; % Actin Retro
    case 4
        noiseParam=[1.96/3.61 0.00031036 2e-4 0.02849789 1.96]; % W512r
    otherwise
        error('noise parameters selection must be 1, 2, 3 or 4');
end

% choose first input directory and its first image
[fileName,dirName] = uigetfile('*.tif','Select the first image of a sequence');
[path,body,no,ext]=getFilenameBody(fileName);

if(isa(fileName,'char') & isa(dirName,'char'))
   % Recover all file names from the stack
   outFileList1=defineStackNames([dirName fileName]);
   % Number of files 
   n=length(outFileList1);
else
   return
end

% adjusting the number of images processed (if needed)
if nbImages>n
    nbImages=n;
end

% extracting the Bit Depth information
aux1=imfinfo([dirName,filesep,fileName]);
BitDepth=14%aux1.BitDepth;
if aux1.ColorType~='grayscale'
    error('Not a Gray Scale Image');
end

% indexing
s=length(num2str(n));
strg=sprintf('%%.%dd',s); 

% choose output directory and file name
[fileOutName,dirOutName] = uiputfile('','Select an output directory and Specify output FILE name down in the Dialog');

% create output sub-directories for the different image groups
mkdir([dirOutName],'NoQuat');
mkdir([dirOutName],'NoTert');
mkdir([dirOutName],'NoSeco');
mkdir([dirOutName],'NoPrim');
mkdir([dirOutName],'Initial');

% processing
h=waitbar(0,'Please wait! The program is substracting the images');
for i=1:nbImages
    
    % create index string
    indxStr=sprintf(strg,i);
    
    if DEBUG==1
        indxStr='01'; %debug
    end
    
    % read input image
    Iinit=imread([dirName,filesep,body,indxStr,ext]);
    
    % initial treatment of the raw data
    I=double(Iinit);
    I=I/(2^BitDepth-1);
    
    % filter image
    IG=gauss2d(I,1);
    
    % set the corner values of the filtered image the same as the corner values of the raw data (low) in order to improve the triangulation
    IG(1,1)=I(1,1);
    IG(1,size(IG,2))=I(1,size(I,2));
    IG(size(IG,1),1)=I(size(I,1),1);
    IG(size(IG,1),size(IG,2))=I(size(I,1),size(I,2));
    
    % local minima
    Imin=locmin2d(IG,[3,3]);
    
    % test primary speckles (intial/filtered) image)
%     [yi,xi,y,x,Imax,candsP]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,enhTriang);  
    [yi,xi,y,x,Imax,candsP,triMin,pMin]=fsmPrepConfirmSpeckles(IG,Imin,noiseParam,enhTriang,userROIbw);
    
    % substracted image
    Inew=fsmPrepSubstructMaxima(IG,Imax,SIG,candsP); 
    
    % test secondary speckles (substracted image)
%     [ysi,xsi,ys,xs,ImaxS,candsS]=fsmPrepConfirmSpeckles(Inew,Imin,noiseParam,enhTriang); 
    [ysi,xsi,ys,xs,ImaxS,candsS]=fsmPrepConfirmLoopSpeckles(Inew,noiseParam,enhTriang,triMin,pMin,IG,userROIbw);
    
    % test for distance criterion
%     [yss,xss,candsS]=fsmPrepCheckDistance(ys,xs,y,x,candsS,candsP,Inew,IG); 
    candsS=fsmPrepCheckDistance(candsS,candsP);
%     ImaxS=fsmPrepUpdateImax(ImaxS,yss,xss);??
    candsTot=cat(2,candsP,candsS);
    
    % twice substracted image
    Inew2=fsmPrepSubstructMaxima(Inew,ImaxS,SIG,candsS); 
    
    % test tertiary speckles (twice substracted image)
%     [yti,xti,yt,xt,ImaxT,candsT]=fsmPrepConfirmSpeckles(Inew2,Imin,noiseParam,enhTriang); 
    [yti,xti,yt,xt,ImaxT,candsT]=fsmPrepConfirmLoopSpeckles(Inew2,noiseParam,enhTriang,triMin,pMin,IG,userROIbw);
    
    % test for distance criterion
    candsT=fsmPrepCheckDistance(candsT,candsTot);
    candsTot=cat(2,candsTot,candsT);
    
    % three times substracted image
    Inew3=fsmPrepSubstructMaxima(Inew2,ImaxT,SIG,candsT); 
    
    % test quatro speckles (three times substracted image)
%     [yqi,xqi,yq,xq,ImaxQ,candsQ]=fsmPrepConfirmSpeckles(Inew3,Imin,noiseParam,enhTriang); 
    [yqi,xqi,yq,xq,ImaxQ,candsQ]=fsmPrepConfirmLoopSpeckles(Inew3,noiseParam,enhTriang,triMin,pMin,IG,userROIbw);
    
    % test for distance criterion
    candsQ=fsmPrepCheckDistance(candsQ,candsTot);
    
    % three times substracted image
    Inew4=fsmPrepSubstructMaxima(Inew3,ImaxQ,SIG,candsQ); 
    
    % cut the borders before saving
    shift=4;  
    Inew4=Inew4(5:end-shift,5:end-5);
    Inew3=Inew3(5:end-shift,5:end-5);
    Inew2=Inew2(5:end-shift,5:end-5);
    Inew=Inew(5:end-shift,5:end-5);
    IG=IG(5:end-shift,5:end-5);
    
    if DEBUG==0
        % convert back the images
        Inew4=Inew4*(2^BitDepth-1);
        Inew3=Inew3*(2^BitDepth-1);
        Inew2=Inew2*(2^BitDepth-1);
        Inew=Inew*(2^BitDepth-1);
        IG=IG*(2^BitDepth-1);
        switch BitDepth
            case 8
                Inew4=uint8(Inew4);
                Inew3=uint8(Inew3);
                Inew2=uint8(Inew2);
                Inew=uint8(Inew);
                IG=uint8(IG);
            otherwise
                Inew4=uint16(Inew4);
                Inew3=uint16(Inew3);
                Inew2=uint16(Inew2);
                Inew=uint16(Inew);
                IG=uint16(IG);
%             otherwise
%                 error('the images are of unknown bit depth!');
        end
        
        % write to file/disk
        imwrite(Inew4,[dirOutName,filesep,'NoQuat',filesep,fileOutName,'NoQuat',indxStr,ext]);
        imwrite(Inew3,[dirOutName,filesep,'NoTert',filesep,fileOutName,'NoTert',indxStr,ext]);
        imwrite(Inew2,[dirOutName,filesep,'NoSeco',filesep,fileOutName,'NoSeco',indxStr,ext]);
        imwrite(Inew,[dirOutName,filesep,'NoPrim',filesep,fileOutName,'NoPrim',indxStr,ext]);
        imwrite(IG,[dirOutName,filesep,'Initial',filesep,fileOutName,'Initial',indxStr,ext]);
    end
    
    waitbar(i/nbImages,h);   
end
close(h);

% debug
if DEBUG==1
    close all
    figure,imshow(IG,[]);
    hold on
    plot(xi-shift,yi-shift,'y*');
    figure,imshow(IG,[]);
    hold on
    plot(x-shift,y-shift,'r.');
    plot(xss-shift,yss-shift,'g.');
    plot(xss-shift,yss-shift,'g.');
    figure,imshow(IG,[]);
    figure,imshow(Inew,[]);
    figure,imshow(Inew2,[]);
    figure,imshow(Inew3,[]);
end