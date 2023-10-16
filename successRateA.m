function SUCCESS_RATE_A = successRateA
% finds the Success Rate of detection of 2GKs next to each other when noise added
% either r or A is fixed plot several curves for SNR 3,5,7 and 8

% close all;
SIG=1.33;
Radi=20;
colordef white;
figure
for i = 1:4
    switch i
        case 1
            SNR = 1e+10;
            colo = 'b';
        case 2
            SNR = 10;
            colo = 'g';
        case 3
            SNR = 7;
            colo = 'r';
        case 4
            SNR = 4;
            colo = 'y';
    end
    % SNR=5; % 3,5,7 or 8!
    r=4.5; % if r fixed - change a (relative intensity) r=4.5 default
    centr = [0 0];
    c=0;
    successR = 0;
    SUCCESS_RATE = [];
    A_plot=[];
    % Initializing progress bar
    h = waitbar(0,'Calculating SUCCESS RATE for the concrete r or ReIn');
    attempts = 2000; % default 1000
    amax = 1;
    amin = 0;
    step = 0.01; % default -0.5
    %--------------------------------------------------------------------------
%     tic;
    % for ReIn=0:1:19%4 %flagPlot=0
    for a=amin:step:amax
        successR=0;
        for AARON = 1:attempts
            I=[];GK1=[];GK2=[];y=[];x=[];
            [j1,GK1]=Gauss2D11center(SIG,Radi,1+a,centr);
            [j1,GK2]=Gauss2D11center(SIG,Radi,1,[centr(1) centr(2)+r]);
            I=zeros(31);
            I=GK1+GK2;    
            I=I+(1/SNR.*randn(2*Radi+1)); % Signal-to-noise ratio
            IG=gauss2d(I,1);
            Imax=locmax2d(IG,[5,5]);      
            [y,x]=find(ne(Imax,0));% find the coordinates/positions of the local maxima      
            c2=0; % counter 
            % the position of the second speckle
            d=createDistanceMatrix([1+Radi+centr(2)+r,1+Radi+centr(1)],[y,x]);
            match = find(d<=1.5);
            if length(match)>=2
                error('two speckles too close')
            end
            if ~isempty(match)
                successR = successR + 1;
            end    
            %             for i=1:length(y)
            %                 if I(y(i),x(i))>0.5% & I(y(i),x(i))<1.1
            %                     c2=c2+1;
            %                 end 
            %             end
            %             if c2 == 2
            %                 successR = successR + 1;
            %             end   
        end
        %Update wait bar
        c=c+1;
        waitbar(c/((amax-amin)/step),h);    
        SUCCESS_RATE=[SUCCESS_RATE,successR/attempts*100];
        A_plot = [A_plot,a+1];
    end
%     toc;
    %Close waitbar
    close(h);
    plot(A_plot,SUCCESS_RATE,colo);
    hold on
end
hold off
ylabel('SUCCESS RATE [%]')
xlabel('RELATIVE INTENSITY')
%--------------------------------------------------------------------------
% DEBUG FIGURES
% figure,imshow(Imax,[]);
% title('I max');
% figure,imshow(IG,[]);
% hold on;
% plot(x,y,'g*');
% plot(1+Radi+centr(1),1+Radi+centr(2)+r,'r.')
% hold off;
% title('the filtered image with the local maxima');




