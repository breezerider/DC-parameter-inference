function ca(varargin)
%% Panel a
%% Comparing results from Cristian's article, Figure 2 in (Mueller, 2012), for CMA-ES and current results using DC

% colloidal aggregation (from Mueller, C. et al., 2012)
tc = 'ca';
k = [2    ...
     2.1  ...
     0.1  ...
     1.0  ...
     0.01 ...
     0.1  ...
     15];
kidx = 2:length(k);

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
inparams.dTimeStart = 3000;
inparams.dTimeStep = 0.1;
inparams.dTimeEnd = 4000;

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

    inparams
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

    % simulation parameters
    ts = inparams.dTimeStart;
    dt = inparams.dTimeStep;
    te = inparams.dTimeEnd;
    [f, t] = costfcnPSSA_CMAES([k(:)' 0 0],tc,ts,dt,te,1); % initialize the cost function
    [f, t] = costfcnPSSA_TTestFano([k(:)' 0 0],tc,ts,dt,te,1,1000,1000); % initialize the cost function
    s = warning('off','MATLAB:load:variableNotFound');
    for k_p = 2 : nargin
      if(2 == exist(varargin{k_p}, 'file'))
        fprintf('Processing %s\n', varargin{k_p});
        load(varargin{k_p},'inparams','outCMA','out','TimeElapsed_CMAES','TimeElapsed_DC');

        k_i = k;
        if(strcmp(inparams.strMethod,'cmaes'))
          k_i(kidx) = 10.^(outCMA.solutions.mean.x); k_i = k_i(:);
          CMAES_avg = sum([CMAES_avg; k_i(kidx)'], 1); CMAES_cnt = CMAES_cnt + 1;
          CMAES_data = [CMAES_data; k_i(kidx)' TimeElapsed_CMAES];
          CMAES_files{end+1} = varargin{k_p};

          fprintf(hFileCMAES, '# %s\n', varargin{k_p});
          for k_pk=1:length(k_i)
            fprintf(hFileCMAES, '%f ', k_i(k_pk));
          end
          f1 = costfcnPSSA_CMAES([k_i(:)' 0 0]);
          [f2, t] = costfcnPSSA_TTestFano([k_i(:)' 0 0]);
          fprintf(hFileCMAES, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
        elseif(strcmp(inparams.strMethod,'dc'))
          k_i(kidx) = 10.^(out.muVec(end, :)); k_i = k_i(:);
          vol = out.volVec(end);
          DC_avg = sum([DC_avg; k_i(kidx)' .* vol], 1); DC_totvol = DC_totvol + vol; DC_avg_plain = sum([DC_avg_plain; k_i(kidx)'], 1); DC_cnt = DC_cnt + 1;
          DC_data = [DC_data; k_i(kidx)' vol TimeElapsed_DC];
          DC_files{end+1} = varargin{k_p};

          fprintf(hFileDC, '# %s, volume of the feasible region = %e\n', varargin{k_p}, vol);
          for k_pk=1:length(k_i)
            fprintf(hFileDC, '%f ', k_i(k_pk));
          end
          f1 = costfcnPSSA_CMAES([k_i(:)' 0 0]);
          [f2, t] = costfcnPSSA_TTestFano([k_i(:)' 0 0]);
          fprintf(hFileDC, '\n# Estimated cost: %f vs (%f, %f)\n', f1, f2, t);
        else
          error('Invalid method')
        end
      else
        fprintf('Warning: file does not exist: ''%s''\n', varargin{k_p});
      end
    end
    warning(s);
    %
    if(~isempty(DC_data))||(~isempty(CMAES_data))
      DC_plot_data = zeros(size(DC_data,1),3); % e, vol, T
      CMAES_plot_data = zeros(size(CMAES_data,1),2); % e, T

      for k_p=1:size(DC_data, 1)
        DC_plot_data(k_p, 1) = sqrt(sum(((DC_data(k_p,1:length(kidx)) - k(kidx)) ./ k(kidx)).^2));
        DC_plot_data(k_p, 2:3) = DC_data(k_p, length(kidx) + [1 2]);
      end

      for k_p=1:size(CMAES_data, 1)
        CMAES_plot_data(k_p, 1) = sqrt(sum(((CMAES_data(k_p,1:length(kidx)) - k(kidx)) ./ k(kidx)).^2));
        CMAES_plot_data(k_p, 2) = CMAES_data(k_p, length(kidx) + 1);
      end

      [~, DC_data_idx] = sort(DC_plot_data, 1);
      [~, CMAES_data_idx] = sort(CMAES_plot_data, 1);

      hFilePLOT_DC = fopen(fullfile(strOutputPath,[strOutputName sprintf('-fig-pe-%s-plot-dc.dat',tc)]),'w+');
      hFilePLOT_CMAES = fopen(fullfile(strOutputPath,[strOutputName sprintf('-fig-pe-%s-plot-cmaes.dat',tc)]),'w+');
      N = 5:5:5*floor(max(size(DC_data,1),size(CMAES_data,1))/5);
      for k_n=1:length(N)
        if(~isempty(DC_data))
          idx = DC_data_idx(1:min(N(k_n),length(DC_data_idx)));
          % correlation (e, Vol)
          R_C1 = corrcoef(DC_plot_data(idx,1), log(DC_plot_data(idx, 2))/log(10));
          %R_C1 = corrcoef(DC_plot_data(idx,1), DC_plot_data(idx, 2));
          R_C1 = R_C1(1,2);
          % correlation (e, T)
          R_C2 = corrcoef(DC_plot_data(idx,1), DC_plot_data(idx, 3));
          R_C2 = R_C2(1,2);
          nval = length(idx) - 2;
          tval = tinv(1.0 - 0.05, nval);
          t1 = R_C1 * sqrt(nval/(1 - R_C1^2));
          t2 = R_C2 * sqrt(nval/(1 - R_C2^2));
          % calculate the regression coeffs for e = f(Vol)
          X = [ones(length(idx), 1) log(DC_plot_data(idx,2))/log(10)];
          %X = [ones(length(idx), 1) DC_plot_data(idx,2)];
          B = X \ DC_plot_data(idx,1);
          %
          fprintf(hFilePLOT_DC, '# Corr coeff DC (n = %i, top %i considered): (e, Vol) = %f, t = %f > %f, p = %f; (e, T) = %f, t = %f > %f, p = %f (t < t_c => R == 0)\n', size(DC_data, 1), length(idx), R_C1, abs(t1), tval, 2*tcdf(-abs(t1), nval), R_C2, abs(t2), tval, 2*tcdf(-abs(t2), nval));
          fprintf(hFilePLOT_DC, '# e = (%f) + (%f) * ln(Vol)/ln(10)\n', B);
          %fprintf(hFilePLOT_DC, '# e = (%f) + (%f) * Vol\n', B);
          clear t1 t2 R_C1 R_C2 tval nval idx
        end
        if(~isempty(CMAES_data))
          idx = CMAES_data_idx(1:min(N(k_n),length(CMAES_data_idx)));
          % correlation (e, T)
          R_C = corrcoef(CMAES_plot_data(idx,1), CMAES_plot_data(idx, 2));
          R_C = R_C(1,2);
          nval = length(idx) - 2;
          tval = tinv(1.0 - 0.05, nval);
          t = R_C * sqrt(nval/(1 - R_C^2));
          tval = tinv(1.0 - 0.05, nval);
          %
          fprintf(hFilePLOT_CMAES, '# Corr coeff CMAES (n = %i, top %i considered): (e, T) = %f, t = %f > %f, p = %f (t < t_c => R == 0)\n', size(CMAES_data, 1), length(idx), R_C, abs(t), tval, 2*tcdf(-abs(t), nval));
          clear t R_C tval nval idx
        end
      end
      clear N k_n
      fmtStr = repmat({'k%i;'},1,length(kidx)); fmtStr = [fmtStr{:}]; fmtStr = fmtStr(1:end-1);
      fprintf(hFilePLOT_DC, ['# ' sprintf(fmtStr,kidx) ';e;Vol;T\n']);
      fprintf(hFilePLOT_CMAES, ['# ' sprintf(fmtStr,kidx) ';e;Vol;T\n']);
      fmtStr = repmat({'%e;'},1,length(kidx)); fmtStr = [fmtStr{:}]; % fmtStr = fmtStr(1:end-1);
      fprintf(hFilePLOT_DC, ['# ground truth : ' sprintf(fmtStr,k(kidx)) '\n']);
      fprintf(hFilePLOT_CMAES, ['# ground truth : ' sprintf(fmtStr,k(kidx)) '\n']);
      for k_p=1:size(DC_plot_data, 1)
        fprintf(hFilePLOT_DC, [ '# ' DC_files{DC_data_idx(k_p)} '\n']);
        fprintf(hFilePLOT_DC, [fmtStr '%e;%e;%e\n'], DC_data(DC_data_idx(k_p), 1:length(kidx)), DC_plot_data(DC_data_idx(k_p), :));
      end
      for k_p=1:size(CMAES_plot_data, 1)
        fprintf(hFilePLOT_CMAES, [ '# ' CMAES_files{CMAES_data_idx(k_p)} '\n']);
        fprintf(hFilePLOT_CMAES, [fmtStr '%e;NaN;%e\n'], CMAES_data(CMAES_data_idx(k_p), 1:length(kidx)), CMAES_plot_data(CMAES_data_idx(k_p), :));
      end
      fclose(hFilePLOT_DC);
      fclose(hFilePLOT_CMAES);
    else
      fprintf('No data\n');
    end

    if isdeployed()
      quit(0)
    else
      return
    end
  end
end

% sample ID
k_s = inparams.nSampleID;

if(strcmp(inparams.strMethod,'cmaes'))
  %% CMAES
  % metric params
  num_moments = 1;

  % produce the target dataset
  clear oraclePSSA_CMAES
  clear costfcnPSSA_TCS_CMAES
  oraclePSSA_CMAES([k(:)' 0 0],'costfcnPSSA_CMAES',kidx,tc,inparams.dTimeStart,inparams.dTimeStep,inparams.dTimeEnd,num_moments);
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
  [t,f] = oraclePSSA_CMA_DC([k(:)' 0 0],'costfcnPSSA_TTestFano',kidx,4 + floor(3 * log(length(kidx))),0.1,tc,inparams.dTimeStart,inparams.dTimeStep,inparams.dTimeEnd,num_moments,num_repetitions,bootstrap_size);

  % DC loop
  %for k_s=1:inparams.nSamples
    % starting point
    xstart = LBounds + rand(length(kidx),1).*(UBounds-LBounds);

    %set up Design Centering options
    clear inopts.initQ
    inopts.LBound        = LBounds;
    inopts.UBound        = UBounds;
    inopts.pn            = 2; %change p-norm
    inopts.MaxEval       = 1000 * length(kidx);
    inopts.SavingModulo  = 100 * length(kidx);
    inopts.VerboseModulo = 100 * length(kidx);
    inopts.hitP_adapt    = 0;

    inopts.Plotting = 'on';
    inopts.unfeasibleSave = 1;

    % oracle options (2)
    inopts.CMA.PopSize = 4 + floor(3 * log(length(kidx)));
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
      inopts.para_hitP_adapt.PVec = 1./exp([1:-0.2:0.2 0.05]);
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
