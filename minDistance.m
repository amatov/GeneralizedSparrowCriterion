function [r_plot,SNR_plot]=minDistance(flagPlot)

Debug=1;
close all;

if nargin==0
    flagPlot=1;
end
sR_plt=100;
SIG=1.33;
Radi=15;
SNR=100000;
ReIn=0;
ReIn_plot=[];
SNR_plot=[];
r_plot=[]; 
c=0;

% Initializing progress bar
% h = waitbar(0,'Calculating critical distances for different SNR');

%for i=1:100

% for ReIn=0:1:19%4 %flagPlot=0
for SNR=15:-1:15% flagPlot=1
    for r=10:-.05:3
        
        [j1,GK1]=Gauss2D11center(SIG,Radi,1+ReIn,[2 2]);
        [j1,GK2]=Gauss2D11center(SIG,Radi,1,[2 2+r]);
        I=zeros(31);
        I=GK1+GK2; 
        
        I=I+(1/SNR.*randn(31)).*flagPlot; % Signal-to-noise ratio(with Flag condition)
        
        I=I./max(I(:));
        I=I.*0.1;%0.038623;% maximum intensity in a real data 
        IG=gauss2d(I,1);
        Imax=locmax2d(IG,[5,5]);
        Imin=locmin2d(IG,[3,3]);
        if flagPlot==1
            % standard deviation
            a=1.4919e-4;
            % mean
            b=0.0282; % SAMO AKO NIaMA NOISE B=0
            % noise parameter
            AB=2e-4;
            %[yi,xi]=find(ne(Imax,0));% all the local maxime before the sifnificance test
            % analyze speckles - validate, locmax, locmin...
            info=fsmPrepBkgEstimDelauNoEnh(size(IG),Imax,Imin);
            [Imax,info]=fsmPrepTestLocalMaxima(IG,Imax,info,[1.96/3.55 a AB b 1.96],IG);
        end
        % find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
        [y,x]=find(ne(Imax,0));
        
        % comparison r and r_hat(d)
%         d=createDistanceMatrix([y(1) x(1)],[y(2) x(2)])
%         r
        
        P=find(ne(Imax,0));
        flagDist=0;% DEBUG
        %P(find(IG(P)<0.5*max(IG(:))))=[]; % TOZI RED GO RAZMARKIRAI SAMO AKO IMA MNOGO NOISE !
        
        if length(P)==1%~isempty(res)
            flagDist=1;% DEBUG
            c=c+1;
            ReIn_plot=[ReIn_plot,ReIn];
            r_plot=[r_plot,r];
            SNR_plot=[SNR_plot,SNR];
            r_break=r% DEBUG
            flagDist;
            break
            
            % Update wait bar
% 	        waitbar(SNR/200,h);
            
        end
         flagDist;% DEBUG
    end
end
%end
SUCCESS_RATE=(100-c)/100
% sR_plt=[sR_plt,SUCCESS_RATE];
% figure,plot(sR_plt);
% Close waitbar
% close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG; FIGURES

if flagPlot==1
    coef=fit_ex(SNR_plot,r_plot) % fitting the results for r(SNR);I1:I2=1:1(fixed)
    
%     figure,plot(SNR_plot,r_plot);
%     xlabel('SNR');
%     ylabel('distance');
    %axis([3 15 3 7]);
    
else
    coef_neg=fit_ex_neg(ReIn_plot,r_plot) % fitting the results for r(I1:I2); SNR=Inf
    
%     figure,plot(ReIn_plot,r_plot);
%     xlabel('I1+" "/I2');
%     ylabel('distance');
%     axis([0 10 3 8]);
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Debug==1
    
    figure,imshow(Imax,[]);
    title('I max');
    figure,imshow(IG,[]);
    hold on;
    plot(x,y,'g*');
    hold off;
    title('the filtered image with the local maxima');
%     figure,imshow(Inew,[]);
%     hold on;
%     plot(xn,yn,'g*');
%     hold off;
%     title('I new - image after substruction');
    figure,surf(I);
    axis([0 31 0 31 -.1 .1]);
    title('unfiltered image');
    figure,surf(IG);
    axis([0 31 0 31 -.1 .1]);
    title('filtered image');
%     figure,surf(IG1);
%     axis([0 31 0 31 -.1 .1]);
%     title('masked with a GK local maxima points with intensity delta.I for the raw data image');
%     figure,surf(Inew);
%     axis([0 31 0 31 -.1 .1]);
%     title('image after substruction');
end
 



