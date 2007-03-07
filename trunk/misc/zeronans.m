%zeronans.m
function y=zeronans(A)
% N=size(A,1);
% y=A;
% for nx=1:N
%     for ny=1:N
%         if isnan(A(nx,ny)) == 1
%             y(nx,ny)=0;
%         end
%     end
% end
y=A;
y(isnan(A))=0;