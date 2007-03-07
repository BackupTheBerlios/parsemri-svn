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

% N = M(M(:,55)~=0,55);
data = find(M(:,55)==8);