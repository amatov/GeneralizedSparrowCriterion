function [r_plot,SNR_plot]=minDistanceSubstruct 
Debug=1;
close all;
SIG=1.88;
Radi=15;
SNR_plot=[];
SNR_substruct_plot=[];
r_plot=[]; 
r_substruct_plot=[];

for SNR=50:-1:5
    for r=(3.74+2.5/SNR):-.05:0.05
        
        [j1,GK1]=Gauss2D11center(SIG,Radi,1,[2 2]);
        [j1,GK2]=Gauss2D11center(SIG,Radi,1,[2 2+r]);
        I=zeros(31);
        I=GK1+GK2; 
        I=I+1/SNR.*randn(31); % NOISE
        I=I./max(I(:));
        I=I.*(0.1-0.0282);%0.038623;% maximum intensity in a real data 
        I=I+0.0282; % dark background
        IG=gauss2d(I,1);
        
%         % set the corner values of the filtered image the same as the corner values of the raw data (low) in order to improve the triangulation
%         IG(1,1)=I(1,1);
%         IG(1,size(IG,2))=I(1,size(I,2));
%         IG(size(IG,1),1)=I(size(I,1),1);
%         IG(size(IG,1),size(IG,2))=I(size(I,1),size(I,2));
         
        d=createDistanceMatrix([2 2],[2 2+r]);% ZASHTO MI TRIABVA ??
        
        Imax=locmax2d(IG,[5,5]);
        Imin=locmin2d(IG,[3,3]);
        % standard deviation
        a=1.4919e-4;
        % mean
        b=0.0282; 
        % noise parameter
        AB=1e-4;
        [yi,xi]=find(ne(Imax,0));% all the local maxime before the sifnificance test
        % analyze speckles - validate, locmax, locmin...
        info=fsmPrepBkgEstimationDelaunay(size(IG),Imax,Imin,1);
        [Imax,info]=fsmPrepTestLocalMaxima(IG,Imax,info,[1.96/3.55 a AB b 1.96],IG);
        
        info;%DEBUG
        
        
        
        % find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
        [y,x]=find(ne(Imax,0));
        P=find(ne(Imax,0));
        flagDist=0;% DEBUG
        P(find(IG(P)<0.5*max(IG(:))))=[];
        
        
        if length(P)==1 % TUI VECHE NE MI TRIABVA !!
            
            flagDist=1;% DEBUG
            r_plot=[r_plot,r];
            SNR_plot=[SNR_plot,SNR];
            r_break=r;% DEBUG
            %break
            
            % THE PART FOR CASE OF ONLY ONE LOCMAX DETECTED - THEN, SUBSTRUCT !
            
            % prepearing Imax for substruction
            Imax=zeros(size(Imax));
            sizey=size(info,3);
            % setting the value of Imax to the values from deltaI, i.e. taking into account the background
            for i=1:sizey
                if info(5,1,i)<0.5*max(IG(:))
                        info(8,1,i)=0;
                end
                if info(8,1,i)==1 & info(6,1,i)>0 % status flag - a local maximun is significant or not          
                    Imax(info(1,1,i),info(1,2,i))=info(6,1,i);
                end

            end
            
            % find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
            [y,x]=find(ne(Imax,0)); % delta Is

            % the masked with a GK local maxima points with intensity delta.I for the raw data image
            IG1=gauss2d1(Imax,SIG);

            % substruction (FROM THE FILTERED IMAGE)
            Inew=IG-IG1;
        
            % new local maxima and minima
            IGnew=Inew;%gauss2d(Inew,1);
            Imaxnew=locmax2d(IGnew,[5,5]);
            Iminnew=locmin2d(IGnew,[3,3]);
 
            % analyze the new speckles - validate, locmax, locmin...
            info=analyzeSpeckles(size(I),Imaxnew,Imin,1);
            [Imaxnew,info]=validateSpeckles2(Inew,Imaxnew,info,[1.96/3.55 a AB b 1.96]);

            % find the coordinates/positions of the local maxima after selecting only the significant local maxima/speckles
            [yn,xn]=find(ne(Imaxnew,0)); 
            
            sizey=size(info,3);
            for i=1:sizey
                if info(5,1,i)<0.9*max(Inew(:))
                        info(8,1,i)=0; % DOBAVI NESHTO KOETO DA GI IZCHISTAVA OT IMAXNEW !!
                        Imaxnew(info(1,1,i),info(1,2,i))=0;
                end
            end
            
            [yn,xn]=find(ne(Imaxnew,0)); 
            
            if length(yn)~=1
                r_substruct_plot=[r_substruct_plot,r];
                SNR_substruct_plot=[SNR_substruct_plot,SNR];
                r_substruct_break=r;
                break
            end
            
        end
          
    end % of the 'r' loop
end % of the SNR loop

coef=fit_ex(SNR_substruct_plot,r_substruct_plot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEBUG

% figure,plot(SNR_plot,r_plot);
% xlabel('SNR');
% ylabel('distance');
% axis([0 50 2 8]);
% 
% figure,plot(SNR_substruct_plot,r_substruct_plot);
% xlabel('SNR');
% ylabel('distance Substruct');
% axis([2 20 0 1.2]);

if Debug==1
    
    % all local minima
    lmy1(1:size(info,3))=info(2,1,1:size(info,3));
    lmx1(1:size(info,3))=info(2,2,1:size(info,3));
    lmy2(1:size(info,3))=info(3,1,1:size(info,3));
    lmx2(1:size(info,3))=info(3,2,1:size(info,3));
    lmy3(1:size(info,3))=info(4,1,1:size(info,3));
    lmx3(1:size(info,3))=info(4,2,1:size(info,3));
    
    figure,imshow(Imax,[]);
    title('I max');
    figure,imshow(IG,[]);
    hold on;
    plot(x,y,'g*');
    plot(xi,yi,'y.');
    hold off;
    title('the filtered image with the local maxima (confirmed ones in green)');
    
    % local maxima and minima + triangulation
    shift=0;
    figure,
    imshow(IG((1+shift):(end-shift),(1+shift):(end-shift)),[]);
    hold on;
    plot(x-shift,y-shift,'r*');
    plot(xi,yi,'y.');
    for v=1:size(lmx1,2)
        plot([lmx1(v)-shift lmx2(v)-shift lmx3(v)-shift lmx1(v)-shift],[lmy1(v)-shift lmy2(v)-shift lmy3(v)-shift lmy1(v)-shift],'w-');
        plot(lmx1-shift,lmy1-shift,'g.');
        plot(lmx2-shift,lmy2-shift,'b.');
        plot(lmx3-shift,lmy3-shift,'y.');
    end
    hold off;
    title('local maxima and minima + triangulation');
    
    figure,imshow(Inew,[]);
    hold on;
    plot(xn,yn,'g*');
    hold off;
    title('I new - image after substruction');
    figure,surf(I);
    axis([0 31 0 31 -.1 .1]);
    title('unfiltered image');
    figure,surf(IG);
    axis([0 31 0 31 -.1 .1]);
    title('filtered image');
    figure,surf(IG1);
    axis([0 31 0 31 -.1 .1]);
    title('masked with a GK local maxima points with intensity delta.I for the raw data image');
    figure,surf(Inew);
    axis([0 31 0 31 -.1 .1]);
    title('image after substruction');
end






