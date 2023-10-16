function SUCCESS_RATE = successRate
% finds the Success Rate of detection of 2GKs next to each other when noise added
% either r or A is fixed plot several curves for SNR 3,5,7 and 8

% close all;
SIG=1.33;
Radi=20;
% colordef black;
figure
for i = 1:4 % default 7:-1:5
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
    % SNR=8; % 3,5,7 or 8!
    ReIn=1; % if ReIn fixed - change r (A=ReIn+1 or A=2)
    centr = [0 0];
    c=0;
    SUCCESS_RATE = [];
    R_plot=[];
    % Initializing progress bar
    h = waitbar(0,'Calculating SUCCESS RATE for the concrete r or ReIn');
    attempts = 2000; % default 1000
    rmax = 5.5;
    rmin = 4;
    step = 0.01; % default -0.05
    %--------------------------------------------------------------------------
%     tic;
    % for ReIn=0:1:19%4 %flagPlot=0
    for r=rmin:step:rmax
        successR=0;
        for AARON = 1:attempts
            I=[];GK1=[];GK2=[];y=[];x=[];
            [j1,GK1]=Gauss2D11center(SIG,Radi,1+ReIn,centr);
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
            %         for i=1:length(y)
            %             if I(y(i),x(i))>0.5 %& I(y(i),x(i))<1.5
            %                 c2=c2+1;
            %             end 
            %         end
            %         if c2 == 2
            %             successR = successR + 1;
            %         end   
        end
        %Update wait bar
        c=c+1;
        waitbar(c/((rmax-rmin)/step),h);    
        SUCCESS_RATE=[SUCCESS_RATE,(successR/attempts)*100];
        R_plot = [R_plot,r];
    end
%     toc;
    %Close waitbar
    close(h);
    plot(R_plot,SUCCESS_RATE,colo);
    hold on
    %-------------------------------------------------
    % cf=polyfit(R_plot,SUCCESS_RATE,6);
    % ffd=cf(1).*R_plot .^6+cf(2).*R_plot.^5+cf(3).*R_plot.^4+cf(4).*R_plot.^3+cf(5).*R_plot.^2+cf(6).*R_plot+cf(7);
    % plot(R_plot,ffd,'r')
    %------------------------------------------
    % x0=[0 0];
    % FIT_LSQ = lsqnonlin(@fit_simp,x0,R_plot,SUCCESS_RATE)
    % plot(R_plot,FIT_LSQ,'r');
    %------------------------------------------------------
end
hold off
ylabel('SUCCESS RATE [%]')
xlabel('DISTANCE')
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




