function [f,t]=oraclePSSA_CMAES(x,cf,xi,varargin)
%%%% INPUT: %%%%
% x     candidate solution OR vector of intial parameters for the model (log10)
% cf    Cost function name
% tc    Test case name
% xi    indexes of the parameters substituted by the candidate solution
% ts    beginning of the target dataset
% dt    time step of the target dataset
% te    final time of the target dataset
% nm    number of moments to consider
%%%% OUTTUT: %%%%
% values of the cost function and the dataset

persistent costfcn kall kidx

if isempty(x)
  error('Parameter vector cannot be empty')
elseif ~isnumeric(x)
  error('Parameter vector must consist of numeric values')
end

if nargin > 2

  % ----- Handle cost function param
  costfcn = cf;
  if isempty(costfcn)
      error('Cost function cannot be empty');
  elseif ~ischar(costfcn) && ~isa(costfcn, 'function_handle')
      error('Cost function argument must be a string or function handle');
  end
  % If given cost function was a string convert it to function handle
  if ~isa(costfcn, 'function_handle')
      costfcn = str2func(costfcn);
  end

  % store all params
  kall = x;

  % param indexes
  if(isempty(xi))
    kidx = 1:length(k);
  else
    kidx = unique(xi(:));
  end

  % a trick to keep it simple below
  x = log10(kall(kidx));
end

ki = kall;
ki(kidx) = 10.^x;

if nargin > 3
  [f, t] = feval(costfcn, ki, varargin{:});
else
  f = feval(costfcn, ki);
end

end
