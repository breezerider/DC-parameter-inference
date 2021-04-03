function [f,T]=compare_TCS_CMAES(x,varargin)
%%%% INTUT: %%%%
% x     candidate solution
% ts    beginning of the target dataset
% dt    time step of the target dataset
% te    final time of the target dataset
% nm    highest order of moments to consider in the cost function
%%%% OUTTUT: %%%%
% values of the cost function and both of its subcomponents

% persistent MU_T ACF_T zx

persistent tstart tstep tend nmoments N n A_POP A_VEC

f = NaN;
T = NaN;

bSave=(nargout > 1);

if nargin > 1
  tstart = varargin{1};
  tstep = varargin{2};
  tend = varargin{3};
  nmoments = varargin{4};
  [N,n,A_POP,A_VEC]=compare_TCS_common;
else
  f1 = 0;
  f2 = 0;

  MU_A = zeros(n, nmoments, 2);
  if bSave
    T = zeros(n, mexpssa_num_timepoints(tend, tstep, tstart), N);
  end

  parfor n_i=1:n
    A = mexpssa('tcs', [x(:)' A_VEC(n_i) A_POP 0 A_POP 0], tend, tstep, tstart);
    if bSave
      T(n_i, :, :) = A;
    end
    MU_A(n_i,:,:) = mexpssa_moments(A(:,[1 5]), 1:nmoments);
  end

  MU_A_den = MU_A(:,:,1);
  MU_A_den(0.0 == MU_A_den) = 1.0;
  f1 = (MU_A(:,:,1) - MU_A(:,:,2)) ./ MU_A_den;

%    fprintf('>>>> Molecule numbers (I -> R): ');
%    R_A_nom = 0;
%    R_A_den = 0;
%    R_A_Aavg= mean(MU_A(:,1));
%    for n_i=1:n
%      MU_Ai = shiftdim(MU_A(n_i,:,:));
%      fprintf('(%i) %f -> %f; ', n_i, MU_Ai(1), MU_Ai(2));
%      R_A_nom = R_A_nom + (MU_Ai(1) - MU_Ai(2)).^2;
%      R_A_den = R_A_den + (R_A_Aavg - MU_Ai(1)).^2;
%    end
%    R_A = 1.0 - R_A_nom/R_A_den;
%    fprintf('R_A = %f\n', R_A);

  f = 2.0 * sum(abs(f1(:)));
end

end
