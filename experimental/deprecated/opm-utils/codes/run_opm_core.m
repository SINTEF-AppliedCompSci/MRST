function output_dir = run_opm_core(param, simulator)
   mrstNargInCheck(0, 2, nargin);

   codes_core = fullfile('examples', 'sim_2p_incomp_reorder');

   nin = nargin;

   if nin < 1, param = []; end

   if nin < 2,
      simulator = simcommand(opm_dir(evalin('base', 'myself')), ...
                             codes_core);
   end

   if ~isempty(param),
      input_filename = fullfile(fileparts(param.deck_filename), ...
                                'test_run.param');
      paramStructToParamFile(param, input_filename);
   end

   command = simulator;
   if ~isempty(param),
      command = [command, ' ', input_filename];
   end

   if mrstVerbose
      status = system(command);
      result = [];
   else
      [status, result] = system(command);
   end

   if status ~= 0,
      disp(result);
      error(['Command failed: ''', command, ''''])

   else
      if isfield(param, 'output_dir'),
         output_dir = param.output_dir;
      else
         output_dir = 'output';
      end
   end
end
