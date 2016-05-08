function clean_output(varargin)

p = inputParser;
p.addParameter('dir', '.');
p.parse(varargin{:});

delete([p.Results.dir '/subdomain*.csv.*']);

end