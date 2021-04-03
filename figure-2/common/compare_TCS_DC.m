function [f,t,T]=compare_TCS_DC(x,varargin)
%%%% INTUT: %%%%
% x     candidate solution
% ts    beginning of the target dataset
% dt    time step of the target dataset
% te    final time of the target dataset

%%%% OUTTUT: %%%%
% values of the cost function and both of its subcomponents

% ts = 4000; dt = 0.1; te = 5000; A = mexpssa('tcs', [ones(6,1)' 1 100 0 100 0], te, dt, ts);
% compare_TCS_DC(ones(6,1),4,0.1,5,1,1000,1000)

persistent tstart tstep tend nmoments nrepetitions bssize N n A_POP A_VEC tval Jval Rc chiVal

bStore=(nargin > 1);
bSave=(nargout > 2);

if bStore
  tstart = varargin{1};
  tstep = varargin{2};
  tend = varargin{3};

  nmoments = varargin{4};
  nrepetitions = varargin{5};
  bssize = varargin{6};

  [N,n,A_POP,A_VEC]=compare_TCS_common;
%    tval = tinv(1.0 - 0.05/(2*nmoments), (n - 1));

  chiVal = chi2inv(1.0 - 0.05, n - 1);

%  ap=ff2n(n); %the all possibilities based on two-level full-factorial design.
%  ap(ap~=1)=-1; %change 0 with -1
%  k=1:1:n;
%  Jval=ap*k'; %all possible sums of ranks for k elements
%  %to compute the p-value see how many values are more extreme of the observed
%  %W and then divide for the total number of combinations

%  tval = tinv(1.0 - 0.05, (n - 2))
%  Rc = tval ./ sqrt(n - 2 + tval.^2)
end

if bSave
  T = zeros(n, mexpssa_num_timepoints(tend, tstep, tstart), N);
end

MU_Aavg = zeros(n, N);
%MU_Astd = zeros(n, N);
parfor n_i=1:n
  A = mexpssa('tcs', [x(:)' A_VEC(n_i) A_POP 0 A_POP 0], tend, tstep, tstart);
  if bSave
%      size(A)
%      A(1, :)
%      A(end, :)
    T(n_i, :, :) = A;
  end

  MU_A = zeros(nrepetitions, nmoments, N);
  for r_i=1:nrepetitions
    A_i = A(randi([1 size(A,1)], bssize, 1),:);
    MU_A(r_i, :, :) = mexpssa_moments(A_i, 1:nmoments);
  end

  MU_Aavg(n_i, :) = shiftdim(mean(MU_A),1);
%    MU_Astd(n_i, :) = shiftdim(std(MU_A),1);
end

%  fprintf('>>>> Molecule numbers (I -> R): ');
%  R_A_nom = 0;
%  R_A_den = 0;
%  R_A_Aavg= mean(MU_Aavg(:,1));
%  for n_i=1:n
%    MU_Ai = shiftdim(MU_Aavg(n_i,[1 5]));
%    fprintf('(%i) %f -> %f; ', n_i, MU_Ai(1), MU_Ai(2));
%    R_A_nom = R_A_nom + (MU_Ai(1) - MU_Ai(2)).^2;
%    R_A_den = R_A_den + (R_A_Aavg - MU_Ai(1)).^2;
%  end
%  R_A = 1.0 - R_A_nom/R_A_den;
%  fprintf('R_A = %f\n', R_A );
%  clear MU_Ai n_i R_A R_A_nom R_A_den

% Linear regression
%  mux = mean(MU_Aavg(:,1));
%  muy = mean(MU_Aavg(:,5));
%  ssq = sum((MU_Aavg(:,1) - mux).^2);
%  b_1 = (sum(MU_Aavg(:,1) .* MU_Aavg(:,5)) - n * mux * muy) / ssq;
%  b_0 = muy - b_1 * mux;
%  fprintf('mux = %f\tmuy = %f\tssq = %f\tb_0 = %f\tb_1 = %f\n', mux, muy, ssq, b_0, b_1);
%  % esq = sum((MU_Aavg(:,1) - MU_Aavg(:,5)).^2)
%  fprintf('>>>> Linear model (I -> R): ');
%  R_A_nom = 0;
%  R_A_den = 0;
%  R_A_Aavg= mean(MU_Aavg(:,1));
%  for n_i=1:n
%    MU_A = shiftdim(MU_Aavg(n_i,[1 5]));
%    fprintf('(%i) %f -> %f; ', n_i, MU_A(1), b_0 + b_1 * MU_A(1));
%    R_A_nom = R_A_nom + (b_0 + b_1 * MU_A(1) - MU_A(2)).^2;
%    R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
%  end
%  R_A = 1.0 - R_A_nom/R_A_den * (n - 1)/(n - 2);
%  fprintf('R_A_adj = %f\n', R_A );
%  clear MU_A n_i R_A R_A_nom R_A_den
%
%  esq = sum((MU_Aavg(:,1) - MU_Aavg(:,5)).^2);
%  sb1 = sqrt((esq / (n - 2)) / ssq );
%  sb0 = sqrt(esq / (n - 2) * (1 / n + mux^2 / ssq));
%  fprintf('esq = %f\t', esq);
%  fprintf('abs(b_0) / sb0 = %f\t', abs(b_0) / sb0);
%  fprintf('abs(b_1 - 1) / sb1 = %f\n', abs(b_1 - 1) / sb1);
%
%  esq = sum(((b_0 + b_1 * MU_Aavg(:,1)) - MU_Aavg(:,5)).^2);
%  sb1 = sqrt((esq / (n - 2)) / ssq );
%  sb0 = sqrt(esq / (n - 2) * (1 / n + mux^2 / ssq));
%  fprintf('esq = %f\t', esq);
%  fprintf('abs(b_0) / sb0 = %f\t', abs(b_0) / sb0);
%  fprintf('abs(b_1 - 1) / sb1 = %f\n', abs(b_1 - 1) / sb1);
%
%  abs(b_0) / sb0;
%  fval = abs(b_1 - 1) / sb1;
%  f = fval;
%  t = tinv(1.0 - 0.05/2.0, (n - 2))

% residual w.r.t. response
%  esq = sum(((b_0 + b_1 * MU_Aavg(:,1)) - MU_Aavg(:,5)).^2);
%  sb1 = sqrt((esq / (n - 2)) / ssq );
%  sb0 = sqrt(esq / (n - 2) * (1 / n + mux^2 / ssq));
%  tval_b0 = abs(b_0) / sb0; % b_0 == 0 ? |t| < t_c => true
%  tval_b1 = abs(b_1) / sb1; % b_1 == 0 ? |t| < t_c => true
% residual w.r.t. inducer
%  esq = sum((MU_Aavg(:,1) - MU_Aavg(:,5)).^2);
%  sb1 = sqrt((esq / (n - 2)) / ssq );
%  tval_b1_I = abs(b_1 - 1) / sb1; % b_1 == 1 ? |t| < t_c => true

%  f = tval_b0;
%  f = max(tval_b0, tval_b1_I);
%  t = tinv(1.0 - 0.01/2.0, (n - 2));
%  fprintf('tval_b0 = %f; t_crit = %f\n', tval_b0, t);
%  fprintf('tval_b0 = %f; tval_b1_I = %f; t_crit = %f\n', tval_b0, tval_b1_I, t);
%  fprintf('tval_b0 = %f; tval_b1 = %f; tval_b1_I = %f; t_crit = %f\n', tval_b0, tval_b1, tval_b1_I, t);

% L-inf norm criterion
%  fval = (MU_Aavg(:,1) - MU_Aavg(:,5)) ./ MU_Aavg(:,1);
%  f = max(abs(fval(:))) * mean(abs(fval(:)))
%  t = max(MU_Astd(:,5) ./ MU_Aavg(:,5))

%  MU_stat = MU_Aavg(:,1) - MU_Aavg(:,5);
% T-test
%  MU_statavg = mean(MU_stat)
%  MU_statstd = std(MU_stat)
%  Tval = MU_statavg ./ MU_statstd * sqrt(n)
%  f = max(abs(Tval(:)))
%  t = tval

% Wilcoxon
%  r=tiedrank(abs(MU_stat)); %ranks and ties
%  W=sum(r.*sign(MU_stat)); %Wilcoxon statics (sum of ranks with sign)
%  p=length(Jval(abs(Jval)>=abs(W)))/length(Jval); %p-value
%  f = p;
%  t = 0.05;

%  R = corrcoef(MU_Aavg(:,1), MU_Aavg(:,5))
%  f = R(1,2) - Rc
%  t = 0

% Chi-Sq
f = sum((MU_Aavg(:,1) - MU_Aavg(:,5)).^2./MU_Aavg(:,1));
t = chiVal;
%  fprintf('chisq: x^2 = %f; x^2_c = %f\n', f, chiVal);

end
