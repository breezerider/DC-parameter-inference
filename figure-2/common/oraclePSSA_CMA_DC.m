function [r,f]=oraclePSSA_CMA_DC(x,cf,xi,n,e,varargin)
%%%% INTUT: %%%%
% x     candidate solution(s)
%       single candidate solution as column vector of log10 of reaction rates
% cf    Cost function name or function handle
% e     threshold for relative error
%%%% OUTTUT: %%%%
% r     feasibility
%   1 if x is feasible
%   0 otherwise
% f     value of the cost function

persistent kall kidx errmarg costfcn threshold past_f bCMA bDC num_evals_per_iter num_evals num_acc num_rej fmtStr

if nargin > 1
  bCMA = true;
  bDC  = false;
num_acc = 0;
num_rej = 0;
  threshold = [];
  past_f = [];

  % ----- Handle cost function param
  if isempty(cf)
    error('Cost function cannot be empty');
  elseif ~ischar(cf) && ~isa(cf, 'function_handle')
    error('Cost function argument must be a string or function handle');
  end
  if ischar(cf) && strcmp(cf,'--reset--')
    fprintf('reset the oracle\n');
    return
  end

  % store all params
  if isempty(x)
    error('Parameter vector cannot be empty')
  elseif ~isnumeric(x)
    error('Parameter vector must consist of numeric values')
  end

  kall = x;

  % param indexes
  if(isempty(xi))
    kidx = 1:length(k);
  else
    kidx = unique(xi(:));
  end

  % If given cost function was a string convert it to function handle
  if ~isa(cf, 'function_handle')
    costfcn = str2func(cf);
  else
    costfcn = cf;
  end

  % number of candidate solutions per iteration of Lp-Adaptation
  if isnumeric(n)
    num_evals_per_iter = round(n);
  else
    error('Invalid value for population size');
  end

  if isnumeric(e)
    errmarg = e;
  else
    errmarg = 0.05
  end

  if nargin > 5
    [f,r] = feval(costfcn, x, varargin{:});
  else
    f = NaN;
    r = NaN;
  end

  fmtStr = repmat({'%e, '},1,length(kidx)); fmtStr = [fmtStr{:}]; fmtStr = fmtStr(1:end-2);

  fprintf('done setting up\n');
  return
end

if isempty(x)
  error('Parameter vector cannot be empty')
elseif ~isnumeric(x)
  error('Parameter vector must consist of numeric values')
end

ki = kall;
ki(kidx) = 10.^x;

[f,t] = feval(costfcn, ki);

if isempty(threshold)||isempty(past_f)
  num_evals = 1;
  threshold = (1 + errmarg/2) * f;
  past_f = zeros(num_evals_per_iter, 1);
  past_f(num_evals) = f;
end

if(false == bDC)
  if(mod(num_evals,num_evals_per_iter)==0)
  %    num_evals = 1;
    if ((1 + errmarg) * mean(past_f(:))) < threshold
      threshold = ((1 + errmarg) * mean(past_f(:)))
    end
  %    past_f = zeros(num_evals_per_iter, 1);
  end

%    t = threshold;

  if(num_evals < num_evals_per_iter)
    num_evals = num_evals + 1;
  else
    num_evals = 1;
  end

  past_f(num_evals) = f;
%  fprintf('Covariance Matrix Adaptation Evolution Strategy\n');
else
  if(num_evals < num_evals_per_iter)
%      t = threshold;
    num_evals = num_evals + 1;
%  fprintf('Covariance Matrix Adaptation Evolution Strategy\n');
  else
%      t = 0;
%      num_evals = num_evals_per_iter + 1;
%      if(true == bCMA)
%        fprintf('Quit Covariance Matrix Adaptation Evolution Strategy, switch to DC\n');
%      end
%      bCMA = false;
%  fprintf('Design centering\n');
  end
  %error('hi')
end

% decision on feasibility
if(f > t*(1 - errmarg)^bCMA)
  if((true == bCMA)&&(f < threshold))
    r = 1;
    fprintf(['CMA: accept point (' fmtStr ') t=%e\n'], ki(kidx(:)), threshold);
  else
    if((false == bCMA)&&(true == bDC))
num_rej = num_rej + 1;
%fmtStr = repmat({'%e, '},1,length(ki)); fmtStr = [fmtStr{:}]; fmtStr = fmtStr(1:end-2);
      fprintf(['DC: reject point #%d (' fmtStr ') f=%f t=%f\n'], num_rej, ki(kidx(:)), f, t);
    elseif(true==bCMA)
      fprintf(['CMA: reject point (' fmtStr ') t=%e\n'], ki(kidx(:)), threshold);
    end
    r = 0;
  end
else
  if(false == bDC)
    fprintf('Prepare for Design Centering\n');
  end
%    else
num_acc = num_acc + 1;
%fmtStr = repmat({'%e, '},1,length(ki)); fmtStr = [fmtStr{:}]; fmtStr = fmtStr(1:end-2);
    fprintf(['DC: accept point #%d (' fmtStr ') f=%f t=%f\n'], num_acc, ki(kidx(:)), f, t);
%    end
  bDC = true;
  bCMA = false;
  r = 1;
end

r = struct('feasible',r,'optimization',bCMA);

end
