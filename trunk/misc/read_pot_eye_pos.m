function [data] = read_plot_eye_pos(FILENAME)
% [data] = read_plot_eye_pos(FILENAME)
%
% Open's eye position data saved in Sam's format and plots
% left/right horizontal/vertical eye position.  returns 4 col mtx with eye
% position

fd = fopen(FILENAME,'r');
M = fread(fd, 'int32',0,'ieee-le');
rows = size(M,1) / 64;
M = reshape(M, 64, rows)';
size(M);

cols = [2 3 5 6 55];
%cols = [2 5 3 6];

titles={'Horizontal Left Eye (positive == right eye mvmt)','Vertical Left Eye (pos == upward eye mvmt)','Horizontal Right Eye (pos == to the right)','Vertical Right Eye (pos == up)'};

%plot(M(:,cols(1)),M(:,cols(2)),M(:,cols(3)),M(:,cols(4)));


%for i=1:4,
%    subplot(2,2,i);
%	hold on
%    plot(M(:,cols(i)));
%    title(titles{i});
%end;

plot(M(:,cols(1)),'Color','r', 'LineStyle','none','Marker','.');
hold on;
plot(M(:,cols(2)),'Color','g', 'LineStyle','none','Marker','.');
hold on;
plot(M(:,cols(3)),'Color','m', 'LineStyle','none','Marker','.');
hold on;
plot(M(:,cols(4)),'Color','b', 'LineStyle','none','Marker','.');
hold on
plot((M(:,55)*1000)-5000,'Color','k');
hold off

Mage=[M(:,cols(1)) M(:,cols(2)) M(:,cols(3)) M(:,cols(4)) (M(:,cols(5))*1000)-5000 ];
%imagesc(Mage');
%image(M');
%surf(M);

data = M(:,cols);
return;
