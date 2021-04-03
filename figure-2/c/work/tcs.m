function tcs(varargin)
%% Panel b
%% Comparing results for the orthodox TCS model (adapted from  Josephines thesis) between CMA-ES (as outlined in (Mueller, 2012)) and current results from DC

% two component bacterial system
tc = 'tcs';
k = ones(1,6);
kidx = 1:length(k);

% bounds on the parameter values
LBounds = floor(log10(k(kidx)) - 1)';
UBounds = floor(log10(k(kidx)) + 1)';
%

if ~isdeployed()
  addpath(genpath('../../common/'));
  addpath(genpath('../../../mex/'));
end

inparams.strOutputPath = ['.'];
inparams.strOutputName = 'Ostrenko2021';
% inparams.nSamples = 10; % number of samples
inparams.nSampleID = 1;
inparams.strMethod = 'dc'; % either 'cmaes' or 'dc'
inparams.dTimeStart = 500;
inparams.dTimeStep = 0.1;
inparams.dTimeEnd = 1000;

if(nargin > 0)
  if(~strcmp(varargin{1},'analyze'))
    valid_inparams = fieldnames(inparams);
    types_inparams = zeros(length(valid_inparams),1);

    for k_p=1:length(valid_inparams)
      if(ischar(inparams.(valid_inparams{k_p})))
        types_inparams(k_p) = 1;
      end
    end

    if mod(nargin,2) ~= 0
      error('Expecting an even number of arguments');
    end
    for k_p = 1 : 2 : nargin-1
      k_s = find(ismember(valid_inparams,varargin{k_p}));
      if isempty(k_s)
        error('Argument "%s" is not a valid parameter name', varargin{k_p});
      end
      if types_inparams(k_s)
        val = varargin{k_p+1};
      else
        val = str2double(varargin{k_p+1});
        if isnan(val)
          error('Value "%s" for variable "%s" is not a numeric scalar', varargin{k_p+1}, varargin{k_p});
        end
      end
      inparams.(varargin{k_p}) = val;
    end
  else
    strOutputPath = inparams.strOutputPath;
    strOutputName = inparams.strOutputName;

    hFileCMAES = fopen(fullfile(strOutputPath,[strOutputName sprintf('-fig-pe-%s-%s.dat','cmaes',tc)]),'w+');
    hFileDC = fopen(fullfile(strOutputPath,[strOutputName sprintf('-fig-pe-%s-%s.dat','dc',tc)]),'w+');
    DC_avg = []; DC_totvol = 0; DC_avg_plain = []; DC_cnt = 0;
    CMAES_avg = []; CMAES_cnt = 0;

    DC_data = [];
    CMAES_data = [];
    DC_files = {};
    CMAES_files = {};

    if nargin < 2
      error('Not enought parameters. MAT-files with simulation results expected.');
    end

    num_workers = str2num(getenv('SLURM_CPUS_ON_NODE'));
    if ~isempty(num_workers)
      parpool(num_workers);
    end
    clear num_workers

    % simulation parameters
    ts = inparams.dTimeStart;
    dt = inparams.dTimeStep;
    te = inparams.dTimeEnd;
    [f, t] = compare_TCS_CMAES(k,ts,dt,te,1); % initialize the cost function
    [f, t] = compare_TCS_DC(k,ts,dt,te,1,1000,1000); % initialize the cost function
    fprintf('\n\nInitialization done...\n\n');
    s = warning('off','MATLAB:load:variableNotFound');
    for k_p = 2 : nargin
      if(2 == exist(varargin{k_p}, 'file'))
        fprintf('Processing %s\n', varargin{k_p});
        load(varargin{k_p},'inparams','outCMA','out','TimeElapsed_CMAES','TimeElapsed_DC');

        if(strcmp(inparams.strMethod,'cmaes'))
          k_i = 10.^(outCMA.solutions.mean.x); k_i = k_i(:);
          CMAES_avg = sum([CMAES_avg; k_i'], 1); CMAES_cnt = CMAES_cnt + 1;
          CMAES_data = [CMAES_data; k_i' TimeElapsed_CMAES];
          CMAES_files{end+1} = varargin{k_p};

          fprintf(hFileCMAES, '# %s\n', varargin{k_p});
          for k_pk=1:length(k_i)
            fprintf(hFileCMAES, '%f ', k_i(k_pk));
          end
          f1 = compare_TCS_CMAES(k_i);
          [f2, t, T] = compare_TCS_DC(k_i);
          fprintf(hFileCMAES, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
          fprintf(hFileCMAES, '# Molecule numbers (I -> R): ');
          R_A_nom = 0;
          R_A_den = 0;
          R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
          for n_i=1:size(T,1)
            MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
            fprintf(hFileCMAES, '(%i) %f -> %f; ', n_i, MU_A(1), MU_A(2));
            R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
            R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
          end
          R_A = 1.0 - R_A_nom/R_A_den;
          fprintf(hFileCMAES, 'R_A = %f\n', R_A);
        elseif(strcmp(inparams.strMethod,'dc'))
          k_i = 10.^(out.muVec(end,:)); k_i = k_i(:);
          vol = out.volVec(end);
          DC_avg = sum([DC_avg; k_i' .* vol], 1); DC_totvol = DC_totvol + vol; DC_avg_plain = sum([DC_avg_plain; k_i'], 1); DC_cnt = DC_cnt + 1;
          DC_data = [DC_data; k_i' vol TimeElapsed_DC];
          DC_files{end+1} = varargin{k_p};

          fprintf(hFileDC, '# %s, volume of the feasible region = %f\n', varargin{k_p}, vol);
          for k_pk=1:length(k_i)
            fprintf(hFileDC, '%f ', k_i(k_pk));
          end
          f1 = compare_TCS_CMAES(k_i);
          [f2, t, T] = compare_TCS_DC(k_i);
          fprintf(hFileDC, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
          fprintf(hFileDC, '# Molecule numbers (I -> R): ');
          R_A_nom = 0;
          R_A_den = 0;
          R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
          for n_i=1:size(T,1)
            MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
            fprintf(hFileDC, '(%i) %f -> %f; ', n_i, MU_A(1), MU_A(2));
            R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
            R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
          end
          R_A = 1.0 - R_A_nom/R_A_den;
          fprintf(hFileDC, 'R_A = %f\n', R_A);
        else
          error('Invalid method')
        end
      else
        fprintf('Warning: file does not exist: ''%s''\n', varargin{k_p});
      end
    end
    warning(s);
    % output the results
    if(~isempty(DC_avg))
      [~, DC_data_idx] = sortrows(DC_data, length(kidx)+1, 'descend');
      DC_avg = DC_avg ./ DC_totvol;
      DC_avg_plain = DC_avg_plain ./ DC_cnt;
    %
      fprintf(hFileDC, '# Plain average\n');
      for k_pk=1:length(DC_avg_plain)
        fprintf(hFileDC, '%f ', DC_avg_plain(k_pk));
      end
      [f1, T] = compare_TCS_CMAES(DC_avg_plain);
      [f2, t, T] = compare_TCS_DC(DC_avg_plain);
      fprintf(hFileDC, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
      fprintf(hFileDC, '# Molecule numbers (I -> R): ');
      R_A_nom = 0;
      R_A_den = 0;
      R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
      for n_i=1:size(T,1)
        MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
        fprintf(hFileDC, '(%i) %f -> %f; ', n_i, MU_A(1), MU_A(2));
        R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
        R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
      end
      R_A = 1.0 - R_A_nom/R_A_den;
      fprintf(hFileDC, 'R_A = %f\n', R_A);
      %
      fprintf(hFileDC, '# Volume-weighted average\n');
      for k_pk=1:length(DC_avg)
        fprintf(hFileDC, '%f ', DC_avg(k_pk));
      end
      [f1, T] = compare_TCS_CMAES(DC_avg);
      [f2, t, T] = compare_TCS_DC(DC_avg);
      fprintf(hFileDC, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
      fprintf(hFileDC, '# Molecule numbers (I -> R): ');
      R_A_nom = 0;
      R_A_den = 0;
      R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
      for n_i=1:size(T,1)
        MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
        fprintf(hFileDC, '(%i) %f -> %f; ', n_i, MU_A(1), MU_A(2));
        R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
        R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
      end
      R_A = 1.0 - R_A_nom/R_A_den;
      fprintf(hFileDC, 'R_A = %f\n', R_A);
      %
      DC_avg = mean(DC_data(DC_data_idx(1:3),1:length(kidx)));
      fprintf(hFileDC, '# Top 3 by volume, average\n');
      for k_pk=1:length(DC_avg)
        fprintf(hFileDC, '%f ', DC_avg(k_pk));
      end
      [f1, T] = compare_TCS_CMAES(DC_avg);
      [f2, t, T] = compare_TCS_DC(DC_avg);
      fprintf(hFileDC, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
      fprintf(hFileDC, '# Molecule numbers (I -> R): ');
      R_A_nom = 0;
      R_A_den = 0;
      R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
      for n_i=1:size(T,1)
        MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
        fprintf(hFileDC, '(%i) %f -> %f; ', n_i, MU_A(1), MU_A(2));
        R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
        R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
      end
      R_A = 1.0 - R_A_nom/R_A_den;
      fprintf(hFileDC, 'R_A = %f\n', R_A);
    end
    if(~isempty(CMAES_avg))
      CMAES_avg = CMAES_avg ./ CMAES_cnt;
      fprintf(hFileCMAES, '# Average\n');
      for k_pk=1:length(CMAES_avg)
        fprintf(hFileCMAES, '%f ', CMAES_avg(k_pk));
      end
      [f1, T] = compare_TCS_CMAES(CMAES_avg);
      [f2, t, T] = compare_TCS_DC(CMAES_avg);
      fprintf(hFileCMAES, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
      fprintf(hFileCMAES, '# Molecule numbers (I -> R): ');
      R_A_nom = 0;
      R_A_den = 0;
      R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
      for n_i=1:size(T,1)
%          size(T(n_i,:,[1 5]))
%          %mean(T(n_i,:,[1 5]))
%          MU_A = shiftdim(mean(shiftdim(T(n_i,:,[1 5]), 1)), 1);
%          size(MU_A)
%          MU_A
        MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
%          size(MU_A)
%          MU_A
        fprintf(hFileCMAES, '(%i) %f -> %f; ', n_i, MU_A(1), MU_A(2));
        R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
        R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
      end
      R_A = 1.0 - R_A_nom/R_A_den;
      fprintf(hFileCMAES, 'R_A = %f\n', R_A);
    end
    fclose(hFileDC);
    fclose(hFileCMAES);

    if(~isempty(DC_data))||(~isempty(CMAES_data))
      DC_plot_data = zeros(size(DC_data,1),3); % R_A, vol, T
      CMAES_plot_data = zeros(size(CMAES_data,1),2); % R_A, T

      for k_p=1:size(DC_data, 1)
        [f2, t, T] = compare_TCS_DC(DC_data(k_p,1:length(kidx)));
        R_A_nom = 0;
        R_A_den = 0;
        R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
        for n_i=1:size(T,1)
          MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
          R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
          R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
        end
        DC_plot_data(k_p, 1) = 1.0 - R_A_nom/R_A_den;
        DC_plot_data(k_p, 2:3) = DC_data(k_p, length(kidx) + [1 2]);
      end

      for k_p=1:size(CMAES_data, 1)
        [f2, t, T] = compare_TCS_DC(CMAES_data(k_p,1:length(kidx)));
        R_A_nom = 0;
        R_A_den = 0;
        R_A_Aavg= (T(1,1,1) + T(size(T,1),1,1)) / 2.0;
        for n_i=1:size(T,1)
          MU_A = shiftdim(mean(T(n_i,:,[1 5])), 1);
          R_A_nom = R_A_nom + (MU_A(1) - MU_A(2)).^2;
          R_A_den = R_A_den + (R_A_Aavg - MU_A(1)).^2;
        end
        CMAES_plot_data(k_p, 1) = 1.0 - R_A_nom/R_A_den;
        CMAES_plot_data(k_p, 2) = CMAES_data(k_p, length(kidx) + 1);
      end

      [~, DC_data_idx] = sort(DC_plot_data, 1, 'descend');
      [~, CMAES_data_idx] = sort(CMAES_plot_data, 1, 'descend');

      hFilePLOT_DC = fopen(fullfile(strOutputPath,[strOutputName sprintf('-fig-pe-%s-plot-dc.dat',tc)]),'w+');
      hFilePLOT_CMAES = fopen(fullfile(strOutputPath,[strOutputName sprintf('-fig-pe-%s-plot-cmaes.dat',tc)]),'w+');
      N = 5:5:5*floor(max(size(DC_data,1),size(CMAES_data,1))/5);
      for k_n=1:length(N)
        if(~isempty(DC_data))
          idx = DC_data_idx(1:min(N(k_n),length(DC_data_idx)));
          % correlation (R_A, Vol)
          R_C1 = corrcoef(DC_plot_data(idx,1), DC_plot_data(idx, 2));
          R_C1 = R_C1(1,2);
          R_C1log = corrcoef(DC_plot_data(idx,1), log(DC_plot_data(idx,2))/log(10));
          R_C1log = R_C1log(1,2);
          % correlation (R_A, T)
          R_C2 = corrcoef(DC_plot_data(idx,1), DC_plot_data(idx, 3));
          R_C2 = R_C2(1,2);
          nval = length(idx) - 2;
          tval = tinv(1.0 - 0.05, nval);
          t1 = R_C1 * sqrt(nval/(1 - R_C1^2));
          t1log = R_C1log * sqrt(nval/(1 - R_C1log^2));
          t2 = R_C2 * sqrt(nval/(1 - R_C2^2));
          % calculate the regression coeffs for R_A = f(Vol)
          X = [ones(length(idx), 1) DC_plot_data(idx,2)];
          Blin = X \ DC_plot_data(idx,1);
          X = [ones(length(idx), 1) log(DC_plot_data(idx,2))/log(10)];
          Blog = X \ DC_plot_data(idx,1);
          %
          fprintf(hFilePLOT_DC, '# Corr coeff DC (n = %i, top %i considered): (R_A, Vol) = %f, t = %f > %f, p = %f; (R_A, ln(Vol)/ln(10)) = %f, t = %f > %f, p = %f; (R_A, T) = %f, t = %f > %f, p = %f (t < t_c => R == 0)\n', size(DC_data, 1), length(idx), R_C1, abs(t1), tval, 2*tcdf(-abs(t1), nval), R_C1log, abs(t1log), tval, 2*tcdf(-abs(t1log), nval), R_C2, abs(t2), tval, 2*tcdf(-abs(t2), nval));
          fprintf(hFilePLOT_DC, '# R_A = (%f) + (%f) * Vol; R_A = (%f) + (%f) * ln(Vol)/ln(10)\n', Blin, Blog);
          clear nval tval t1 t2 R_C1 R_C2 X Blin Blog
        end
        if(~isempty(CMAES_data))
          idx = CMAES_data_idx(1:min(N(k_n),length(CMAES_data_idx)));
          % correlation (e, T)
          R_C = corrcoef(CMAES_plot_data(idx,1), CMAES_plot_data(idx, 2));
          R_C = R_C(1,2);
          nval = length(idx) - 2;
          tval = tinv(1.0 - 0.05, nval);
          t = R_C * sqrt(nval/(1 - R_C^2));
          %
          fprintf(hFilePLOT_CMAES, '# Corr coeff CMAES (n = %i, top %i considered): (R_A, T) = %f, t = %f > %f, p = %f (t < t_c => R == 0)\n', size(CMAES_data, 1), length(idx), R_C, abs(t), tval, 2*tcdf(-abs(t), nval));
          clear nval tval t R_C
        end
      end
      fmtStr = repmat({'k%i;'},1,length(kidx)); fmtStr = [fmtStr{:}]; % fmtStr = fmtStr(1:end-1);
      fprintf(hFilePLOT_DC, ['# ' sprintf(fmtStr,kidx) 'R_A;Vol;T\n']);
      fprintf(hFilePLOT_CMAES, ['# ' sprintf(fmtStr,kidx) 'R_A;Vol;T\n']);
      fmtStr = repmat({'%e;'},1,length(kidx)); fmtStr = [fmtStr{:}]; % fmtStr = fmtStr(1:end-1);
      for k_p=1:size(DC_plot_data, 1)
        fprintf(hFilePLOT_DC, [ '# ' DC_files{DC_data_idx(k_p)} '\n']);
        fprintf(hFilePLOT_DC, [fmtStr '%e;%e;%e\n'], DC_data(DC_data_idx(k_p), 1:length(kidx)), DC_plot_data(DC_data_idx(k_p), :));
      end
      for k_p=1:size(CMAES_plot_data, 1)
        fprintf(hFilePLOT_CMAES, [ '# ' CMAES_files{CMAES_data_idx(k_p)} '\n']);
        fprintf(hFilePLOT_CMAES, [fmtStr '%e;NaN;%e\n'], CMAES_data(CMAES_data_idx(k_p), 1:length(kidx)), CMAES_plot_data(CMAES_data_idx(k_p), :));
      end
%        fprintf(hFilePLOT, '# R_A;Vol;T\n');
%        for k_p=1:size(DC_plot_data, 1)
%          fprintf(hFilePLOT, '%e;%e;%e\n', DC_plot_data(k_p, :));
%        end
%        for k_p=1:size(CMAES_plot_data, 1)
%          fprintf(hFilePLOT, '%e;NaN;%e\n', CMAES_plot_data(k_p, :));
%        end
      fclose(hFilePLOT_DC); fclose(hFilePLOT_CMAES);
    end

    if isdeployed()
      quit(0)
    else
      return
    end
  end
end

inparams

% sample ID
k_s = inparams.nSampleID;

num_workers = str2num(getenv('SLURM_CPUS_ON_NODE'));
if ~isempty(num_workers)
  parpool(num_workers);
end
clear num_workers

%  aprfor i=1:str2num(getenv('SLURM_CPUS_ON_NODE'))
%    disp('test')
%  end

if(strcmp(inparams.strMethod,'cmaes'))
  %% CMAES
  % metric params
  num_moments = 1;

  % produce the target dataset
  clear oraclePSSA_CMAES
  clear costfcnPSSA_TCS_CMAES
  oraclePSSA_CMAES(k,'compare_TCS_CMAES',kidx,inparams.dTimeStart,inparams.dTimeStep,inparams.dTimeEnd,num_moments);
  %
  % CMA-ES loop
  %for k_s=1:inparams.nSamples
    % starting point
    xstart = LBounds + rand(length(kidx),1).*(UBounds-LBounds);

    % Initial sigma
    initSigma = 1;

    % set up CMAES options
    clear inopts
    inopts.MaxFunEvals = 1000 * length(k);
    inopts.StopFitness = 1e-9;
    inopts.IncPopSize = 2;
    inopts.Restarts = 10;
    inopts.LBounds = LBounds;
    inopts.UBounds = UBounds;
  %    inopts.Display = 'on';
    inopts.Noise.on = 1;
  %    inopts.LogPlot = 'on';

    TimeElapsed_CMAES = tic
    [xmin,fmin,counteval,stopflag,outCMA,bestever] = cmaes_355('oraclePSSA_CMAES',xstart,initSigma,inopts);
    TimeElapsed_CMAES = toc(TimeElapsed_CMAES)

    save(fullfile(inparams.strOutputPath,[inparams.strOutputName sprintf('-pe-cmaes-%s-sample-%i.mat',tc,k_s)]));
  %end
elseif(strcmp(inparams.strMethod,'dc'))
  %% DC
  % metric params
  num_moments = 1;
  num_repetitions = 1000;
  bootstrap_size = 1000;

  % oracle options (1)
  oracle='oraclePSSA_CMA_DC';
  [t,f] = oraclePSSA_CMA_DC(k,'compare_TCS_DC',kidx,4 + floor(3 * log(length(k))),0.1,inparams.dTimeStart,inparams.dTimeStep,inparams.dTimeEnd,num_moments,num_repetitions,bootstrap_size);

%  figure; hold on; Names=cell(2*2,1); for n_i=1:2 for s_i=[1 5] plot(t(n_i,:,s_i)); end; Names((n_i-1)*2 + [1 2]) = {sprintf('population %i - indicator',n_i); sprintf('population %i - response',n_i)}; end; hold off; legend(Names);
%  error('hi')

  % DC loop
  %for k_s=1:inparams.nSamples
    % starting point
    xstart = LBounds + rand(length(kidx),1).*(UBounds-LBounds);

    %set up Design Centering options
    clear inopts.initQ
    inopts.LBound        = LBounds;
    inopts.UBound        = UBounds;
    inopts.pn            = 2; %change p-norm
    inopts.MaxEval       = 1000 * length(k);
    inopts.SavingModulo  = 100 * length(k);
    inopts.VerboseModulo = 100 * length(k);
    inopts.hitP_adapt    = 0;

    inopts.Plotting = 'on';
    inopts.unfeasibleSave = 1;

    % oracle options (2)
    inopts.CMA.PopSize = 4 + floor(3 * log(length(k)));
    %inopts.oracleInopts = {ts,dt,te,num_moments,num_repetitions,bootstrap_size,e,inopts.CMA.PopSize};
    inopts.nOut = 2;
    inopts.valP = 0.5;
    inopts.hitP_adapt = 0;

    oraclePSSA_CMA_DC([],'--reset--');

    TimeElapsed_DC = tic
    out = LpAdaptation(oracle,xstart,inopts);

    if isfield(out,'doDC') && out.doDC
      fprintf('Doing DC!\n');
      inopts.initQ = out.initQ;
      inopts.hitP_adapt = 1;
      inopts.para_hitP_adapt.PVec = 1./exp([1:-0.1:0.5]);
      inopts.para_hitP_adapt.maxEvalSchedule = [1/2 repmat(1/10, 1, 5)];%
      inopts.para_hitP_adapt.numLastSchedule = [1/2 repmat(3/4, 1, 5)];

      out = LpAdaptation(oracle,out.xstart,inopts);
    else
      fprintf('Converged in a single run!\n');
    end
    TimeElapsed_DC = toc(TimeElapsed_DC)

    save(fullfile(inparams.strOutputPath,[inparams.strOutputName sprintf('-pe-dc-%s-sample-%i.mat',tc,k_s)]));
  %end
else
  error('Invalid method')
end

if isdeployed()
  quit(0)
end

end
