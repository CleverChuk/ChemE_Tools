function [antCoff1, antCoff2] = getAntione( name, name2, varargin )
%Gets antoine coefficient for two components


antCoff = load('AntoineTable.mat');
antCoff = antCoff.AntoineTable;

compName = antCoff(:,1);
found = 0;
%{
    TODO
    Implement more than two look ups
%}

if nargin > 2
    %FLAG = remove once implemented
    return;
end
for i = 1:length(compName)
    
    if strcmpi(compName(i), name)
        found = 1;
        antCoff1 = antCoff(i,(2:4));
        antCoff1 = cell2mat(antCoff1);      
    end
end
if ~found
    disp('Component 1 not found')
end
    
for i = 1:length(compName)

    if strcmpi(compName(i), name2)
        found = 1;
        antCoff2 = antCoff(i,(2:4));
        antCoff2 = cell2mat(antCoff2);     
    end
end
if ~found
      disp('Component 2 not found')
end

end

