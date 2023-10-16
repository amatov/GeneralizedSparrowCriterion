function LifeTest


speckleArray(1).status='b';
speckleArray(2).status='s';
speckleArray(3).status='d';
speckleArray(4).status='b';
speckleArray(5).status='s';
speckleArray(6).status='d';
speckleArray(7).status='b';
speckleArray(8).status='s';
speckleArray(9).status='s';
speckleArray(10).status='d';


ll=length(speckleArray)-1
% Mean life time of a speckle
i=1; LifeTime=0; Position=0; count=0; meanLifeTime=[];
for i=1:ll
	if speckleArray(i).status=='b'
        i=i+1;
        count=1;
		while speckleArray(i).status~='d' & speckleArray(i+1).status~='b'
            if i==length(speckleArray)
                break
            end
            i=i+1;
            count=count+1;	
%             if speckleArray(i).status=='b'
%                 i=i-1;
%                 count=0;
%                 break
%             end         
        end
        Position=Position+1;
        LifeTime(Position)=count;
        count=0;
    end
%     i=i+1;
end
meanLifeTime=mean(LifeTime)
