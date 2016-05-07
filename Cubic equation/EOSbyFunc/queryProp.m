clc; clear all;
query = input('Compound name: ','s');
props = load ('props.mat'); % Loads property data
props = props.props;
compounds = props(:,2); % extracts compound names

for i = 1:length(compounds)
            if strcmpi(compounds(i),query)
                rowNum = i
                break
            end
end