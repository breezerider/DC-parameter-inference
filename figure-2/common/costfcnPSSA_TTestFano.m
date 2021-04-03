function [f,t]=costfcnPSSA_TTestFano(x,varargin)
%%%% INTUT: %%%%
% x     candidate solution(s)
%       single candidate solution as column vector of log10 of reaction rates
%%%% OUTTUT: %%%%
% f     T-test statistic value
% t     T-test threshold

persistent testcase tstart tstep tend nmoments nrepetitions bssize MU_T MU_Tavg bFano N

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

  if nargin < 5
    nmoments = 1;
  elseif isnumeric(varargin{5})
    nmoments = round(varargin{5});
  else
    error('Number of moments must be a numeric value');
  end

  if nargin < 6
    nrepetitions = 1;
  elseif isnumeric(varargin{6})
    nrepetitions = round(varargin{6});
  else
    error('Number of repetitions must be a numeric value');
  end

  if nargin < 7
    bssize = 0.1 * mexpssa_num_timepoints(tend, tstep, tstart);
  elseif isnumeric(varargin{7})
    bssize = round(varargin{7});
  else
    error('Bootstrap size must be a numeric value');
  end

  bFano = strcmp(testcase,'ed') || strcmp(testcase,'ca') || strcmp(testcase,'homo') || strcmp(testcase,'sbd'); % 

  if(strcmp(testcase,'homo') || strcmp(testcase,'sbd'))
    N = 1;
  elseif(strcmp(testcase,'ca') || strcmp(testcase,'clc'))
    N = x(1);
  elseif(strcmp(testcase,'ed'))
    N = 4;
  elseif(strcmp(testcase,'tcs'))
    N = 5;
  end
end

A = mexpssa(testcase, x, tend, tstep, tstart);

if(1 == nrepetitions)
  MU_A = mexpssa_moments(A, 1:nmoments);
else
  MU_A = zeros(nrepetitions, nmoments, N);
  for r_i=1:nrepetitions
    A_i = A(randi([1 size(A,1)], bssize, 1),:);
    MU_A(r_i, :, :) = mexpssa_moments(A_i, 1:nmoments);
  end
end

if nargin < 2

  % compute the Fano factors
  %
  Fano_Var = [];
  if bFano
    disp('Fano correction factor used to estimate variance');
    % compute the species population at steady state
    problem.options = optimoptions(@lsqnonlin,'FiniteDifferenceType','central','SpecifyObjectiveGradient',true,'Display','none'); % ,'Algorithm','levenberg-marquardt'
    problem.solver = 'lsqnonlin';
    problem.x0 = MU_Tavg(1,:);
    problem.lb = zeros(N, 1);
    problem.objective = @(s) mexpssa_odes(testcase,s,x);

    [S_ss,resnorm,residual,exitflag,output] = lsqnonlin(problem);

    if(any((S_ss < 1).*(S_ss > 0)))
      S_ss
      disp('Invalid value of the discrete steady state: fractioned molecules detected!');
      %bFano = 0
    else
  %    if(abs(S_ss - MU_Tavg(1,:))./MU_Tavg(1,:) > 0.25)
  %      S_ss = MU_Tavg(1,:);
  %    end

      [J, Q] = mexpssa_lyapunov(testcase, S_ss, x)
      C = lyap(J, Q);

  %    fprintf('<mu> = %f; ss = %f; <dlt_mu^2> = %f\n', MU_Tavg(1,:), S_ss, diag(C));

      % Variance as estimated by LNA
      Fano_Var = sqrt(MU_Tavg(1,:) ./ S_ss .* diag(C)');
    end
  end

  MU_stat = MU_T - MU_A;

  %  fprintf('E[MU_stat]\n');
  if(1 == nrepetitions)
    MU_statavg = mean(MU_stat);
  else
    MU_statavg = shiftdim(mean(MU_stat),1);
  end
  %  fprintf('Var[MU_stat]\n');
  if(1 == nrepetitions)
    MU_statstd = std(MU_stat);
  else
    MU_statstd = shiftdim(std(MU_stat),1);
  end

  if ~isempty(Fano_Var)
    MU_statstd(1,:) = Fano_Var;
  end

  T = MU_statavg ./ MU_statstd;

  f = max(abs(T(:)));

  t = tinv(1.0 - 0.05/(2*nmoments*N), (nrepetitions - 1));

else
  MU_T = MU_A;
  %  fprintf('E[MU_A] before\n');
  MU_Tavg = shiftdim(mean(MU_T),1)
  %  fprintf('Var[MU_A] before\n');
  %MU_Tstd = shiftdim(std(MU_T),1);
end

end
