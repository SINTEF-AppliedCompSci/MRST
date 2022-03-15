classdef mrstSettingsGUI < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                matlab.ui.Figure
        


        TextArea_OutputDir      matlab.ui.control.TextArea
        BrowseButton_OutputDir  matlab.ui.control.Button
        
        TextArea_DataDir        matlab.ui.control.TextArea
        BrowseButton_DataDir    matlab.ui.control.Button
        
        TextArea_allowDL        matlab.ui.control.TextArea
        TextArea_allowDLDesc    matlab.ui.control.TextArea
        ButtonGroup_allowDL     matlab.ui.container.ButtonGroup
        Button_OnAllowDL        matlab.ui.control.RadioButton
        Button_OffAllowDL       matlab.ui.control.RadioButton        
        
        TextArea_promptDL       matlab.ui.control.TextArea
        TextArea_promptDLDesc   matlab.ui.control.TextArea
        ButtonGroup_promptDL    matlab.ui.container.ButtonGroup
        Button_OnPromptDL       matlab.ui.control.RadioButton
        Button_OffPromptDL      matlab.ui.control.RadioButton
        
        TextArea_promptMEX      matlab.ui.control.TextArea
        TextArea_promptMEXDesc  matlab.ui.control.TextArea        
        ButtonGroup_promptMEX   matlab.ui.container.ButtonGroup
        Button_OnPromptMEX      matlab.ui.control.RadioButton
        Button_OffPromptMEX     matlab.ui.control.RadioButton

        TextArea_useMEX         matlab.ui.control.TextArea
        TextArea_useMEXDesc     matlab.ui.control.TextArea
        ButtonGroup_useMEX      matlab.ui.container.ButtonGroup
        Button_OnUseMEX         matlab.ui.control.RadioButton
        Button_OffUseMEX        matlab.ui.control.RadioButton

        TextArea_useOMPDesc     matlab.ui.control.TextArea        
        TextArea_useOMP         matlab.ui.control.TextArea
        ButtonGroup_useOMP      matlab.ui.container.ButtonGroup
        Button_OnUseOMP         matlab.ui.control.RadioButton
        Button_OffUseOMP        matlab.ui.control.RadioButton
        
        TextArea_useHashDesc    matlab.ui.control.TextArea        
        TextArea_useHash        matlab.ui.control.TextArea
        ButtonGroup_useHash     matlab.ui.container.ButtonGroup
        Button_OnUseHash        matlab.ui.control.RadioButton
        Button_OffUseHash       matlab.ui.control.RadioButton        

        OKButton                matlab.ui.control.Button
        CancelButton            matlab.ui.control.Button
        ApplyButton             matlab.ui.control.Button
        HelpButton              matlab.ui.control.Button
        
        ClickheretochangemrstDataDirectoryLabel  matlab.ui.control.Label
        ClickheretochangemrstOutputDirectoryLabel  matlab.ui.control.Label
    end

    
    properties (Access = public)
        settingsHandle % Description
        tempSettings % Store settings temporarily
    end
    

    methods (Access = private)
% 
        % Code that executes after component creation
        function startupFcn(app)
            app.UIFigure.Visible = 'on';
        end

        % Button pushed function: BrowseButton_OutputDir
        function BrowseButton_OutputDirPushed(app, event)
            selpath = uigetdir(app.tempSettings.outputDirectory,'Choose mrstOutputDirectory');
            if selpath == 0
                msg = {'This is the directory where MRST will save generated output such as simulation results.'; ...
                ''; 'Directory has not been changed. mrstOutputDirectory: '; app.tempSettings.outputDirectory};
                
            else
                app.tempSettings.outputDirectory = selpath;                
                msg = {'This is the directory where MRST will save generated output such as simulation results.'; ...
                ''; 'mrstOutputDirectory: '; app.tempSettings.outputDirectory};         
            end
            
            
            app.TextArea_OutputDir.Value = msg; 
            app.UIFigure.Visible = 'off';
            app.UIFigure.Visible = 'on';            
        
        end

        % Button pushed function: BrowseButton_DataDir
        function BrowseButton_DataDirPushed(app, event)
            selpath = uigetdir(app.tempSettings.dataDirectory,'Choose mrstDataDirectory');
            if selpath == 0
                msg = {'This is the directory where MRST will save downloaded datasets (geological models, grids, etc) for simulations.';...
                ''; 'Directory has not been changed. mrstDataDirectory: '; app.tempSettings.dataDirectory};
            else
                app.tempSettings.dataDirectory = selpath;                
                msg = {'This is the directory where MRST will save downloaded datasets (geological models, grids, etc) for simulations.';...
                ''; 'mrstDataDirectory: '; app.tempSettings.dataDirectory};         
            end


            app.TextArea_DataDir.Value = msg; 
            app.UIFigure.Visible = 'off';
            app.UIFigure.Visible = 'on';    

        end

        % Button pushed function: ApplyButton
        function ApplyButtonPushed(app, event)
            updateSettingsStruct(app);
        end

        % Button pushed function: CancelButton
        function CancelButtonPushed(app, event)
            delete(app);
        end

        % Button pushed function: OKButton
        function OKButtonPushed(app, event)
            updateSettingsStruct(app);
            delete(app);
        end

        % Button pushed function: HelpButton
        function HelpButtonPushed(app, event)
            message = {'MRST has several advanced configuration settings.';
                'If you do not want to set these up now, reasonable defaults will be used.';
                ' Settings can be configured later by calling ''mrstSettings'' from the command prompt.'};
            uialert(app.UIFigure,message,'mrstSettings Help','Icon','info');
        end

        % Selection changed function: ButtonGroup_promptDL
        function ButtonGroup_promptDLSelectionChanged(app, event)
            switch app.ButtonGroup_promptDL.SelectedObject
                case app.Button_OnPromptDL
                    app.tempSettings.promptDL.value = true;
                    app.tempSettings.promptDL.defaulted = true;
                case app.Button_OffPromptDL
                    app.tempSettings.promptDL.value = false;
                    app.tempSettings.promptDL.defaulted = false;
            end
        end

        % Selection changed function: ButtonGroup_allowDL
        function ButtonGroup_allowDLSelectionChanged(app, event)
            switch app.ButtonGroup_allowDL.SelectedObject
                case app.Button_OnAllowDL
                    app.tempSettings.allowDL.value = true;
                    app.tempSettings.allowDL.defaulted = true;
                case app.Button_OffAllowDL
                    app.tempSettings.allowDL.value = false;
                    app.tempSettings.allowDL.defaulted = false;
            end
        end

        % Selection changed function: ButtonGroup_promptMEX
        function ButtonGroup_promptMEXSelectionChanged(app, event)
            switch app.ButtonGroup_promptMEX.SelectedObject
                case app.Button_OnPromptMEX
                    app.tempSettings.promptMEX.value = true;
                    app.tempSettings.promptMEX.defaulted = false;
                case app.Button_OffPromptMEX
                    app.tempSettings.promptMEX.value = false;
                    app.tempSettings.promptMEX.defaulted = true;
            end            
        end

        % Selection changed function: ButtonGroup_useMEX
        function ButtonGroup_useMEXSelectionChanged(app, event)
            switch app.ButtonGroup_useMEX.SelectedObject
                case app.Button_OnUseMEX
                    app.tempSettings.useMEX.value = true;
                    app.tempSettings.useMEX.defaulted = true;
                case app.Button_OffUseMEX
                    app.tempSettings.useMEX.value = false;
                    app.tempSettings.useMEX.defaulted = false;
            end            
        end

        % Selection changed function: ButtonGroup_useOMP
        function ButtonGroup_useOMPSelectionChanged(app, event)
            switch app.ButtonGroup_useOMP.SelectedObject
                case app.Button_OnUseOMP
                    app.tempSettings.useOMP.value = true;
                    app.tempSettings.useOMP.defaulted = true;
                case app.Button_OffUseOMP
                    app.tempSettings.useOMP.value = false;
                    app.tempSettings.useOMP.defaulted = false;
            end            
        end

        % Selection changed function: ButtonGroup_useHash
        function ButtonGroup_useHashSelectionChanged(app, event)
            switch app.ButtonGroup_useHash.SelectedObject
                case app.Button_OnUseHash
                    app.tempSettings.useHash.value = true;
                    app.tempSettings.useHash.defaulted = false;
                case app.Button_OffUseHash
                    app.tempSettings.useHash.value = false;
                    app.tempSettings.useHash.defaulted = true;
            end            
        end        
        function updateSettingsStruct(app)
            app.settingsHandle.outputDirectory = app.tempSettings.outputDirectory;
            app.settingsHandle.dataDirectory = app.tempSettings.dataDirectory;
            app.settingsHandle.promptDL = app.tempSettings.promptDL;
            app.settingsHandle.allowDL = app.tempSettings.allowDL;
            app.settingsHandle.promptMEX = app.tempSettings.promptMEX;
            app.settingsHandle.useMEX = app.tempSettings.useMEX;
            app.settingsHandle.useOMP = app.tempSettings.useOMP;  
            app.settingsHandle.useHash = app.tempSettings.useHash;             
        end
        
    end

    % App initialization and construction
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure
            app.UIFigure = uifigure('visible','off');
            app.UIFigure.Position = [100 100 879 859];
            app.UIFigure.Name = 'MRST Settings';
            
            % Set directories
            
            % For easily adjusting positions
            newY = 129;
            
            % Create TextArea_OutputDir
            app.TextArea_OutputDir = uitextarea(app.UIFigure);
            app.TextArea_OutputDir.Editable = 'off';
            app.TextArea_OutputDir.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_OutputDir.Position = [21 641+newY 840 70];
            app.TextArea_OutputDir.Value = {'This is the directory where MRST will save generated output such as simulation results.'; ...
                ''; 'Current mrstOutputDirectory: '; app.settingsHandle.outputDirectory};

            % Create ClickheretochangemrstOutputDirectoryLabel
            app.ClickheretochangemrstOutputDirectoryLabel = uilabel(app.UIFigure);
            app.ClickheretochangemrstOutputDirectoryLabel.Position = [531 609+newY 228 22];
            app.ClickheretochangemrstOutputDirectoryLabel.Text = 'Click here to change mrstOutputDirectory';                        
            
            % Create BrowseButton_OutputDir
            app.BrowseButton_OutputDir = uibutton(app.UIFigure, 'push');
            app.BrowseButton_OutputDir.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_OutputDirPushed, true);
            app.BrowseButton_OutputDir.Position = [761 609+newY 100 22];
            app.BrowseButton_OutputDir.Text = 'Browse';
            
            % Create TextArea_DataDir
            app.TextArea_DataDir = uitextarea(app.UIFigure);
            app.TextArea_DataDir.Editable = 'off';
            app.TextArea_DataDir.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_DataDir.Position = [21 521+newY 840 70];
              app.TextArea_DataDir.Value = {'This is the directory where MRST will save downloaded datasets (geological models, grids, etc) for simulations.';...
                ''; 'Current mrstDataDirectory: '; app.settingsHandle.dataDirectory};
            
            % Create ClickheretochangemrstDataDirectoryLabel
            app.ClickheretochangemrstDataDirectoryLabel = uilabel(app.UIFigure);
            app.ClickheretochangemrstDataDirectoryLabel.Position = [531 489+newY 220 22];
            app.ClickheretochangemrstDataDirectoryLabel.Text = 'Click here to change mrstDataDirectory';            

            % Create BrowseButton_DataDir
            app.BrowseButton_DataDir = uibutton(app.UIFigure, 'push');
            app.BrowseButton_DataDir.ButtonPushedFcn = createCallbackFcn(app, @BrowseButton_DataDirPushed, true);
            app.BrowseButton_DataDir.Position = [761 489+newY 100 22];
            app.BrowseButton_DataDir.Text = 'Browse';      
            
            



            
            % allowDL
            
            % Create TextArea_allowDL
            app.TextArea_allowDL = uitextarea(app.UIFigure);
            app.TextArea_allowDL.Editable = 'off';
            app.TextArea_allowDL.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_allowDL.Position = [21 441+newY 150 30];
            app.TextArea_allowDL.Value = {'Allow file download'};

            % Create ButtonGroup_allowDL
            app.ButtonGroup_allowDL = uibuttongroup(app.UIFigure);
            app.ButtonGroup_allowDL.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_allowDLSelectionChanged, true);
            app.ButtonGroup_allowDL.Position = [181 411+newY 130 60];

            % Create Button_OnAllowDL
            app.Button_OnAllowDL = uiradiobutton(app.ButtonGroup_allowDL);
            app.Button_OnAllowDL.Text = 'On (default)';
            app.Button_OnAllowDL.Position = [11 37 85 22];
            app.Button_OnAllowDL.Value = app.tempSettings.allowDL.value;

            % Create Button_OffAllowDL
            app.Button_OffAllowDL = uiradiobutton(app.ButtonGroup_allowDL);
            app.Button_OffAllowDL.Text = 'Off';
            app.Button_OffAllowDL.Position = [11 7 38 22];
            if app.tempSettings.allowDL.value
                app.Button_OffAllowDL.Value = false;
            else 
                app.Button_OffAllowDL.Value = true;
            end                            

            % Create TextArea_allowDLDesc
            app.TextArea_allowDLDesc = uitextarea(app.UIFigure);
            app.TextArea_allowDLDesc.Editable = 'off';
            app.TextArea_allowDLDesc.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_allowDLDesc.Position = [321 391+newY 540 80];
            app.TextArea_allowDLDesc.Value = {'Turn this on if you want MRST to download data sets and other utilities from the internet. '; ''; 'You can still manually activate downloads through optional inputs to specific functions even when this is switched off.'; ''};

            
            
            % promptDL
            
            % Create TextArea_promptDL
            app.TextArea_promptDL = uitextarea(app.UIFigure);
            app.TextArea_promptDL.Editable = 'off';
            app.TextArea_promptDL.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_promptDL.Position = [21 341+newY 150 40];
            app.TextArea_promptDL.Value = {'Prompt before downloading files'};            
            
            % Create ButtonGroup_promptDL
            app.ButtonGroup_promptDL = uibuttongroup(app.UIFigure);
            app.ButtonGroup_promptDL.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_promptDLSelectionChanged, true);
            app.ButtonGroup_promptDL.Position = [181 321+newY 130 60];

            % Create Button_OnPromptDL
            app.Button_OnPromptDL = uiradiobutton(app.ButtonGroup_promptDL);
            app.Button_OnPromptDL.Text = 'On (default)';
            app.Button_OnPromptDL.Position = [11 37 85 22];

            % Create Button_OffPromptDL
            app.Button_OffPromptDL = uiradiobutton(app.ButtonGroup_promptDL);
            app.Button_OffPromptDL.Text = 'Off';
            app.Button_OffPromptDL.Position = [11 7 38 22];
            if app.tempSettings.promptMEX.value
                app.Button_OffPromptDL.Value = false;
            else 
                app.Button_OffPromptDL.Value = true;
            end                        

            % Create TextArea_promptDLDesc
            app.TextArea_promptDLDesc = uitextarea(app.UIFigure);
            app.TextArea_promptDLDesc.Editable = 'off';
            app.TextArea_promptDLDesc.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_promptDLDesc.Position = [321 321+newY 540 60];
            app.TextArea_promptDLDesc.Value = {'Turn this on if you want to be prompted before MRST downloads any files.'; ''; ''};

            
            % useMEX
            
            % Create TextArea_useMEXDesc
            app.TextArea_useMEXDesc = uitextarea(app.UIFigure);
            app.TextArea_useMEXDesc.Editable = 'off';
            app.TextArea_useMEXDesc.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_useMEXDesc.Position = [321 201+newY 540 110];
            app.TextArea_useMEXDesc.Value = {'Turn this on if you want to use MEX-extensions to improve performance.'; ''; 'Using compiled MEX-extensions can speed up parts of MRST. This requires a C/C++ compiler to be set up with Matlab (see ''''mex -setup'''' for more details). '; ''; 'When this is turned on MRST will attempt to compile MEX files when needed. Building of such files can still be triggered manually.'};

            % Create ButtonGroup_useMEX
            app.ButtonGroup_useMEX = uibuttongroup(app.UIFigure);
            app.ButtonGroup_useMEX.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_useMEXSelectionChanged, true);
            app.ButtonGroup_useMEX.Position = [181 251+newY 130 60];

            % Create Button_OnUseMEX
            app.Button_OnUseMEX = uiradiobutton(app.ButtonGroup_useMEX);
            app.Button_OnUseMEX.Text = 'On (default)';
            app.Button_OnUseMEX.Position = [11 37 85 22];

            % Create Button_OffUseMEX
            app.Button_OffUseMEX = uiradiobutton(app.ButtonGroup_useMEX);
            app.Button_OffUseMEX.Text = 'Off';
            app.Button_OffUseMEX.Position = [11 7 38 22];
            if app.tempSettings.useMEX.value
                app.Button_OffUseMEX.Value = false;
            else 
                app.Button_OffUseMEX.Value = true;
            end                  

            % Create TextArea_useMEX
            app.TextArea_useMEX = uitextarea(app.UIFigure);
            app.TextArea_useMEX.Editable = 'off';
            app.TextArea_useMEX.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_useMEX.Position = [21 271+newY 150 40];
            app.TextArea_useMEX.Value = {'Use MEX-extensions'};
            

            % promptMEX
            
            % Create TextArea_promptMEX
            app.TextArea_promptMEX = uitextarea(app.UIFigure);
            app.TextArea_promptMEX.Editable = 'off';
            app.TextArea_promptMEX.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_promptMEX.Position = [21 151+newY 150 40];
            app.TextArea_promptMEX.Value = {'Prompt before building MEX-extensions'};
            
            % Create ButtonGroup_promptMEX
            app.ButtonGroup_promptMEX = uibuttongroup(app.UIFigure);
            app.ButtonGroup_promptMEX.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_promptMEXSelectionChanged, true);
            app.ButtonGroup_promptMEX.Position = [181 131+newY 130 60];

            % Create Button_OnPromptMEX
            app.Button_OnPromptMEX = uiradiobutton(app.ButtonGroup_promptMEX);
            app.Button_OnPromptMEX.Text = 'On';
            app.Button_OnPromptMEX.Position = [11 37 38 22];

            % Create Button_OffPromptMEX
            app.Button_OffPromptMEX = uiradiobutton(app.ButtonGroup_promptMEX);
            app.Button_OffPromptMEX.Text = 'Off (default)';
            app.Button_OffPromptMEX.Position = [11 7 85 22];
            if app.tempSettings.promptMEX.value
                app.Button_OffPromptMEX.Value = false;
            else 
                app.Button_OffPromptMEX.Value = true;
            end            

            % Create TextArea_promptMEXDesc
            app.TextArea_promptMEXDesc = uitextarea(app.UIFigure);
            app.TextArea_promptMEXDesc.Editable = 'off';
            app.TextArea_promptMEXDesc.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_promptMEXDesc.Position = [321 131+newY 540 60];
            app.TextArea_promptMEXDesc.Value = {'Turn this on if you want to be prompted before MRST builds any MEX-extensions.'; ''; ''};



            % useOMP            
            
            % Create TextArea_useOMP
            app.TextArea_useOMP = uitextarea(app.UIFigure);
            app.TextArea_useOMP.Editable = 'off';
            app.TextArea_useOMP.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_useOMP.Position = [21 81+newY 150 40];
            app.TextArea_useOMP.Value = {'Use OMP when using MEX acceleration'};

            % Create ButtonGroup_useOMP
            app.ButtonGroup_useOMP = uibuttongroup(app.UIFigure);
            app.ButtonGroup_useOMP.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_useOMPSelectionChanged, true);
            app.ButtonGroup_useOMP.Position = [181 61+newY 130 60];

            % Create Button_OnUseOMP
            app.Button_OnUseOMP = uiradiobutton(app.ButtonGroup_useOMP);
            app.Button_OnUseOMP.Text = 'On (default)';
            app.Button_OnUseOMP.Position = [11 37 85 22];

            % Create Button_OffUseOMP
            app.Button_OffUseOMP = uiradiobutton(app.ButtonGroup_useOMP);
            app.Button_OffUseOMP.Text = 'Off';
            app.Button_OffUseOMP.Position = [11 7 38 22];
            if app.tempSettings.useOMP.value
                app.Button_OffUseOMP.Value = false;
            else 
                app.Button_OffUseOMP.Value = true;
            end            

            % Create TextArea_useOMPDesc
            app.TextArea_useOMPDesc = uitextarea(app.UIFigure);
            app.TextArea_useOMPDesc.Editable = 'off';
            app.TextArea_useOMPDesc.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_useOMPDesc.Position = [321 51+newY 540 70];
            app.TextArea_useOMPDesc.Value = {'Turn this on if you want to use OpenMP when using MEX-extensions. '; ''; 'Many MEX extensions include OpenMP support for parallelization. Turn this on if you have OpenMP installed.'};

            
            % useHash            
            
            % Create TextArea_useHash
            app.TextArea_useHash = uitextarea(app.UIFigure);
            app.TextArea_useHash.Editable = 'off';
            app.TextArea_useHash.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_useHash.Position = [21 130 150 40];
            app.TextArea_useHash.Value = {'Use hashing'};            
            
            % Create ButtonGroup_useHash
            app.ButtonGroup_useHash = uibuttongroup(app.UIFigure);
            app.ButtonGroup_useHash.SelectionChangedFcn = createCallbackFcn(app, @ButtonGroup_useHashSelectionChanged, true);            
            app.ButtonGroup_useHash.Position = [181 110 130 60];

            % Create Button_OnUseHash
            app.Button_OnUseHash = uiradiobutton(app.ButtonGroup_useHash);
            app.Button_OnUseHash.Text = 'On';
            app.Button_OnUseHash.Position = [11 37 38 22];

            % Create Button_OffUseHash
            app.Button_OffUseHash = uiradiobutton(app.ButtonGroup_useHash);
            app.Button_OffUseHash.Text = 'Off (default)';
            app.Button_OffUseHash.Position = [11 7 85 22];
            if app.tempSettings.useHash.value
                app.Button_OffUseHash.Value = false;
            else 
                app.Button_OffUseHash.Value = true;
            end                  

            % Create TextArea_useOMPDesc_2
            app.TextArea_useHashDesc = uitextarea(app.UIFigure);
            app.TextArea_useHashDesc.Editable = 'off';
            app.TextArea_useHashDesc.BackgroundColor = [0.9412 0.9412 0.9412];
            app.TextArea_useHashDesc.Position = [321 50 540 120];
            app.TextArea_useHashDesc.Value = {'Turn this on if you want to use hashing for consistency checks and comparing simulation setups.'; ''; 'MRST will use md5 checksums of simulation setups in the AD-OO framework. This makes it easy to detect differences between setups, verify that discrete operators are consistent with the grid and rock, etc. '; ''; 'This feature is largely experimental at present and may increase overheads when running many models. The average user is advised to leave this off.'};

            
            % End buttons
            
            % Create OKButton
            app.OKButton = uibutton(app.UIFigure, 'push');
            app.OKButton.ButtonPushedFcn = createCallbackFcn(app, @OKButtonPushed, true);
            app.OKButton.Position = [511 19 80 22];
            app.OKButton.Text = 'OK';

            % Create CancelButton
            app.CancelButton = uibutton(app.UIFigure, 'push');
            app.CancelButton.ButtonPushedFcn = createCallbackFcn(app, @CancelButtonPushed, true);
            app.CancelButton.Position = [601 19 80 22];
            app.CancelButton.Text = 'Cancel';

            % Create ApplyButton
            app.ApplyButton = uibutton(app.UIFigure, 'push');
            app.ApplyButton.ButtonPushedFcn = createCallbackFcn(app, @ApplyButtonPushed, true);
            app.ApplyButton.Position = [691 19 80 22];
            app.ApplyButton.Text = 'Apply';

            % Create HelpButton
            app.HelpButton = uibutton(app.UIFigure, 'push');
            app.HelpButton.ButtonPushedFcn = createCallbackFcn(app, @HelpButtonPushed, true);
            app.HelpButton.Position = [781 19 80 22];
            app.HelpButton.Text = 'Help';


        end
    end

    methods (Access = public)

        % Construct app
        function app = mrstSettingsGUI(settings)
            
            app.settingsHandle = settings;

            
            app.tempSettings.outputDirectory = app.settingsHandle.outputDirectory;
            app.tempSettings.dataDirectory = app.settingsHandle.dataDirectory;
            app.tempSettings.promptDL = app.settingsHandle.promptDL;
            app.tempSettings.allowDL = app.settingsHandle.allowDL;
            app.tempSettings.promptMEX = app.settingsHandle.promptMEX;
            app.tempSettings.useMEX = app.settingsHandle.useMEX;
            app.tempSettings.useOMP = app.settingsHandle.useOMP;
            app.tempSettings.useHash = app.settingsHandle.useHash;
            
            % Create and configure components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @(app)startupFcn(app))

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
