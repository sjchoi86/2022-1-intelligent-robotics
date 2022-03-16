function str = vec2str(vec,format)
% vector to string
%
% [1 2 3] => [1, 2, 3]
%
if nargin == 1
    format = '%.4f';
end
str = strjoin(cellstr(num2str(vec(:),format)),', ');
str = ['[',str,']'];
