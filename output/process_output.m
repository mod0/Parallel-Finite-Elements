function stepData = process_output(varargin)

p = inputParser;
p.addParameter('dir', '.');
p.addParameter('subdomain', -1, @(x) x >= -1);
p.parse(varargin{:});

subdomain = p.Results.subdomain;
if subdomain == -1
    subdomain = '[0-9]+';
else
    subdomain = num2str(subdomain);
end

maxTokenVal = @(t) max(cell2mat(cellfun(@(x) str2double(x{:}), vertcat(t{:}), 'UniformOutput', false)));

files = dir(p.Results.dir);
tokens = regexp({files.name}, ['^subdomain' subdomain '.csv.([0-9]+)$'], 'tokens');

numSteps = maxTokenVal(tokens) + 1;

stepData = cell(numSteps, 1);

for i = 1:numSteps
    subdomainFiles = regexp({files.name}, ['^subdomain' subdomain '.csv.' num2str(i - 1) '$'], 'match');
    subdomainData = cell2mat(cellfun(@csvread, vertcat(subdomainFiles{:}), 'UniformOutput', false));
    stepData{i} = sortrows(subdomainData, [1 2]);
end

end