close all
FilesS1 = dir('C:\tudatBundle\tudatApplications\Thesis\WP1\PropagationResults\WP1_Segment1_*.dat');

Nr = size(FilesS1);

t=10000;

for( i = 1:Nr )
    clear x y
    Tempdata = importdata(FilesS1(i).name);
    x = Tempdata(:,20);
    y = Tempdata(:,21);
    %plot(x(1:t),y(1:t))
    plot(x,y,'color','b')
    hold on
end

FilesS2 = dir('C:\tudatBundle\tudatApplications\Thesis\WP1\PropagationResults\WP1_Segment2_*.dat');

clear Nr;
Nr = size(FilesS2);

t=10000;

for( j = 1:Nr )
    j;
    clear x y
    Tempdata = importdata(FilesS2(j).name);
    x = Tempdata(:,20);
    y = Tempdata(:,21);
    %plot(x(1:t),y(1:t))
    plot(x,y,'color','r')
    hold on
end

FilesS3 = dir('C:\tudatBundle\tudatApplications\Thesis\WP1\PropagationResults\WP1_Segment3_*.dat');

clear Nr;
Nr = size(FilesS3);

t=100;

for( k = 1:Nr )
    k;
    clear x y
    Tempdata = importdata(FilesS3(k).name);
    x = Tempdata(:,20);
    y = Tempdata(:,21);
    %plot(x(1:t),y(1:t))
    plot(x,y,'color','g')
    hold on
end

MoonX = Tempdata(:,14);
MoonY = Tempdata(:,15);
plot(MoonX,MoonY,'color','k')

plot(0,0,'*')
axis([-8E8 18E8 -10E8 10E8])