function [f,t]=costfcnPSSA_CMAES(x,varargin)
%%%% INTUT: %%%%
% x     candidate solution(s)
%       single candidate solution as column vector of log10 of reaction rates
% tc    Test case name, must be either 'clc' or 'ca'
% xall  vector of reaction rates for the chemical reaction network
% xidx  indexes of the reaction rates to be substituted by the candidate solution
% T     target result
%       single trajectory of produced using sought-for reactions
% ts    beginning of the target dataset
% dt    time step of the target dataset
% te    final time of the target dataset
% MU_T  moments of the target dataset
% ACF_T auto-correlation function
%%%% OUTTUT: %%%%
% values of the cost function and both of its subcomponents

persistent testcase tstart tstep tend nmoments MU_T ACF_T zx

f = NaN;
t = NaN;

if nargin > 1
  if isempty(varargin{1}) || ~ischar(varargin{1})
    error('Test case must be a non-empty string id');
  else
    testcase = varargin{1};
  end

  if isnumeric(varargin{2})
    tstart = varargin{2};
  else
    error('Start time must be a numeric value');
  end

  if isnumeric(varargin{3})
    tstep = varargin{3};
  else
    error('Time step must be a numeric value');
  end

  if isnumeric(varargin{4})
    tend = varargin{4};
  else
    error('End time must be a numeric value');
  end

  if isnumeric(varargin{5})
    nmoments = round(varargin{5});
  else
    error('Number of moments must be a numeric value');
  end
end

A = mexpssa(testcase, x, tend, tstep, tstart);

if nargin > 1
  MU_T = mexpssa_moments(A, 1:nmoments);

  zx = 0;
  ACF_T = zeros(size(A,1)-1, size(A,2));
  ZCV = zeros(1, size(A,2));
  for l=0:(size(A,1)-1)
  %    fprintf('lag %i\n', l);
    ACF_T(l+1,:) = mexpssa_acf(A,l);

    if(0 == l)
      ZCV = sign(ACF_T(l+1,:));
    elseif(0 == zx)
      if(any(ZCV - sign(ACF_T(l+1,:))))
        zx = l;
        break
      end
    else
      break
    end
  end

  ACF_T = ACF_T(1:(zx+1), :);
else
  MU_A = mexpssa_moments(A, 1:nmoments);

  MU_T_den = MU_T;
  MU_T_den(0.0 == MU_T_den) = 1.0;
  f1 = (MU_A - MU_T) ./ MU_T_den;

  ACF_A = mexpssa_acf(A,0:size(ACF_T,1)-1);

  ACF_NUM = sum((ACF_T - ACF_A).^2);
  ACF_DEN = sum((ACF_T + ACF_A).^2);

  f2 = sqrt(ACF_NUM ./ ACF_DEN);

  f = 2 * (sum(sqrt(sum(f1.^2, 2))) + sum(f2));
end

end
