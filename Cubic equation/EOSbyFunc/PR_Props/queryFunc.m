function [ num ] = queryFunc( query )
% Gets compound row number
%
props = load ('props.mat'); % Loads property data
props = props.props;
compounds = props(:,2); % extracts compound names

for i = 1:length(compounds)
    if strcmpi(compounds(i),query)
        num = i;
        break
    end
end


end

