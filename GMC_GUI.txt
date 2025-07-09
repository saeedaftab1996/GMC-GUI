function Geomechanics_GUI
    % Geomechanics_GUI creates a MATLAB GUI that imports geophysical log data from an Excel
    % file, computes stress components (total and effective using Biot's theory), brittleness,
    % elastic (Lame) parameters, and displays four tabs with results, advanced analysis,
    % advanced statistics (including uncertainty), and sensitivity analysis.
    %
    % Save this file as GDA_1_GUI.m and run it in MATLAB (R2018a or later).

    %% Create UI figure and TabGroup
    fig = uifigure('Name', 'Geophysical Stress & Brittleness Calculator', ...
                   'Position', [50 50 1000 900]);
    
    tg = uitabgroup(fig, 'Position', [10 10 980 880]);

    % ------------------------------
    % Tab 1: Main Results
    tabMain = uitab(tg, 'Title', 'Main Results');

    btnLoad = uibutton(tabMain, 'push', ...
        'Text', 'Load Excel Data', ...
        'Position', [20, 820, 150, 30], ...
        'ButtonPushedFcn', @(btn, event) loadDataCallback());

    btnReplot = uibutton(tabMain, 'push', ...
        'Text', 'Replot Results', ...
        'Position', [20, 780, 150, 30], ...
        'ButtonPushedFcn', @(btn, event) updateStressPlot());

    pnlPlotOpts = uipanel(tabMain, 'Title', 'Stress Plot Options', ...
                          'Position', [20, 600, 150, 160]);
    chkVertical = uicheckbox(pnlPlotOpts, 'Text', 'Vertical Stress', ...
                             'Position', [10, 110, 130, 20], 'Value', true);
    chkShMin = uicheckbox(pnlPlotOpts, 'Text', 'Min Horizontal Stress', ...
                          'Position', [10, 85, 130, 20], 'Value', true);
    chkShMax = uicheckbox(pnlPlotOpts, 'Text', 'Max Horizontal Stress', ...
                          'Position', [10, 60, 130, 20], 'Value', true);
    chkPore = uicheckbox(pnlPlotOpts, 'Text', 'Pore Pressure', ...
                         'Position', [10, 35, 130, 20], 'Value', true);
    chkEff = uicheckbox(pnlPlotOpts, 'Text', 'Effective Stresses', ...
                        'Position', [10, 10, 130, 20], 'Value', true);

    resultsTable = uitable(tabMain, 'Position', [200, 600, 760, 240], ...
        'ColumnName', {'Depth (m)', 'Sv (MPa)', 'Pp (MPa)', 'Sh_min (MPa)', ...
                       'Sh_max (MPa)', 'Diff. (MPa)', 'Brittleness (%)'});

    axStress = uiaxes(tabMain, 'Position', [200, 50, 550, 530]);
    title(axStress, 'Stress vs Depth');
    axEff = uiaxes(tabMain, 'Position', [760, 50, 200, 530]);
    title(axEff, 'Effective Stress vs Depth');
    axBrittle = uiaxes(tabMain, 'Position', [20, 50, 160, 530]);
    title(axBrittle, 'Brittleness vs Depth');

    % ------------------------------
    % Tab 2: Advanced Analysis with Scrollable Panel
    tabAdvanced = uitab(tg, 'Title', 'Advanced Analysis');
    btnReplotAdvanced = uibutton(tabAdvanced, 'push', ...
        'Text', 'Replot Advanced', ...
        'Position', [20, 860, 150, 30], ...
        'ButtonPushedFcn', @(btn, event) updateAdvancedPlots());

    advScrollPanel = uipanel(tabAdvanced, 'Position', [20, 20, 940, 820]);
    advScrollPanel.Scrollable = 'on';
    advContentPanel = uipanel(advScrollPanel, 'Position', [0, 0, 940, 1100]);

    axLame = uiaxes(advContentPanel, 'Position', [20, 600, 450, 400]);
    title(axLame, 'Scatter Plot: \lambda vs \mu');
    axScatterUCS = uiaxes(advContentPanel, 'Position', [500, 600, 450, 400]);
    title(axScatterUCS, 'UCS vs Effective Stress');
    axScatterHoriz = uiaxes(advContentPanel, 'Position', [20, 100, 450, 400]);
    title(axScatterHoriz, 'Horizontal Effective Stress vs Strain');
    axScatterVert = uiaxes(advContentPanel, 'Position', [500, 100, 450, 400]);
    title(axScatterVert, 'Scatter Plot: Tensile Strength vs Porosity');

    % ------------------------------
    % Tab 3: Advanced Statistics
    tabStats = uitab(tg, 'Title', 'Advanced Statistics');
    btnUpdateStats = uibutton(tabStats, 'push', ...
        'Text', 'Update Statistics', ...
        'Position', [20, 820, 150, 30], ...
        'ButtonPushedFcn', @(btn, event) updateStatsTable());
    statsTable = uitable(tabStats, 'Position', [20, 450, 940, 300]);
    axErrorStats = uiaxes(tabStats, 'Position', [20, 50, 940, 380]);
    title(axErrorStats, 'Error Plot: Mean \pm Std Dev');

    % ------------------------------
    % Tab 4: Sensitivity Analysis
    tabSensitivity = uitab(tg, 'Title', 'Sensitivity Analysis');
    btnUpdateSensitivity = uibutton(tabSensitivity, 'push', ...
        'Text', 'Update Sensitivity Analysis', ...
        'Position', [20, 820, 200, 30], ...
        'ButtonPushedFcn', @(btn, event) updateSensitivityAnalysis());
    axSensitivity = uiaxes(tabSensitivity, 'Position', [20, 450, 940, 300]);
    title(axSensitivity, 'Sensitivity Analysis (Coefficient of Variation)');
    sensitivityTable = uitable(tabSensitivity, 'Position', [20, 50, 940, 380]);

    %% Store computed data in appdata (initially empty)
    setappdata(fig, 'ComputedData', []);

    %% ----- NESTED CALLBACK FUNCTIONS -----

    function loadDataCallback()
        [file, path] = uigetfile('*.xlsx', 'Select Geophysical Data Excel File');
        if isequal(file, 0)
            uialert(fig, 'No file selected.', 'File Selection');
            return;
        end
        filename = fullfile(path, file);
        try
            data = readtable(filename);
        catch ME
            uialert(fig, ['Error reading Excel file: ', ME.message], 'File Error');
            return;
        end

        if ~any(strcmpi(data.Properties.VariableNames, 'Depth'))
            n = height(data);
            data.Depth = (1:n)';
        end

        reqFields = {'Vp','Vs','Density','Porosity','GammaRay','Resistivity'};
        for i = 1:length(reqFields)
            if ~any(strcmpi(data.Properties.VariableNames, reqFields{i}))
                uialert(fig, ['Missing required field: ', reqFields{i}], 'Data Error');
                return;
            end
        end

        depth       = data.Depth;
        Vp          = data.Vp;
        Vs          = data.Vs;
        Density     = data.Density;
        Porosity    = data.Porosity;
        GammaRay    = data.GammaRay;
        Resistivity = data.Resistivity;
        if any(strcmpi(data.Properties.VariableNames, 'ShaleVolume'))
            ShaleVolume = data.ShaleVolume;
        else
            ShaleVolume = zeros(height(data),1);
        end

        %% Constants & Biot's theory parameter
        g = 9.81;
        water_density = 1000;
        dz = [0; diff(depth)];

        % Total (lithostatic) vertical stress (Pa)
        S_v = (2000*g .* 2000) + cumsum(Density*1000 .* g .* dz);

        % Hydrostatic pore pressure (Pa)
        Pp = water_density * g .* depth;

        % Biot coefficient (default = 1.0; adjust as needed)
        alpha = 1.0;

        %% Poisson's Ratio
        ratio2 = (Vp ./ Vs).^2;
        nu = (ratio2 - 2) ./ (2*(ratio2 - 1));
        nu(isnan(nu)|isinf(nu)) = 0.25;

        %% Horizontal stresses
        S_h_min = S_v .* (nu./(1-nu));
        delta_sigma = 0.3 * S_v;
        S_h_max = S_h_min + delta_sigma;
        diff_stress = S_h_max - S_h_min;

        %% Effective stresses (Biot)
        eff_S_v     = S_v       - alpha*Pp;
        eff_S_h_min = S_h_min   - alpha*Pp;
        eff_S_h_max = S_h_max   - alpha*Pp;
        eff_diff    = eff_S_h_max - eff_S_h_min;

        %% Brittleness
        E_app = Density*1000 .* (Vp.^2);
        E_norm = (E_app - min(E_app)) ./ (max(E_app) - min(E_app) + eps);
        GR_norm = (GammaRay - min(GammaRay)) ./ (max(GammaRay) - min(GammaRay) + eps);
        brittleness = 100 * E_norm .* (1 - GR_norm);

        %% Lame parameters
        mu     = Density*1000 .* (Vs.^2);
        lambda = Density*1000 .* (Vp.^2) - 2*mu;

        %% Convert to MPa
        S_v_MPa        = S_v       /1e6;
        Pp_MPa         = Pp        /1e6;
        S_h_min_MPa    = S_h_min   /1e6;
        S_h_max_MPa    = S_h_max   /1e6;
        diff_MPa       = diff_stress/1e6;
        eff_S_v_MPa    = eff_S_v   /1e6;
        eff_S_h_min_MPa= eff_S_h_min/1e6;
        eff_S_h_max_MPa= eff_S_h_max/1e6;
        eff_diff_MPa   = eff_diff  /1e6;

        %% Store in struct
        computedData = struct( ...
            'depth', depth, ...
            'S_v_MPa', S_v_MPa, ...
            'Pp_MPa', Pp_MPa, ...
            'S_h_min_MPa', S_h_min_MPa, ...
            'S_h_max_MPa', S_h_max_MPa, ...
            'diff_MPa', diff_MPa, ...
            'eff_S_v_MPa', eff_S_v_MPa, ...
            'eff_S_h_min_MPa', eff_S_h_min_MPa, ...
            'eff_S_h_max_MPa', eff_S_h_max_MPa, ...
            'eff_diff_MPa', eff_diff_MPa, ...
            'brittleness', brittleness, ...
            'ShaleVolume', ShaleVolume, ...
            'E_app', E_app, ...
            'mu', mu, ...
            'lambda', lambda, ...
            'Porosity', Porosity, ...
            'BiotCoeff', alpha);
        setappdata(fig, 'ComputedData', computedData);

        %% Update table & plots
        results = table(depth, S_v_MPa, Pp_MPa, S_h_min_MPa, S_h_max_MPa, ...
                        diff_MPa, brittleness, ...
            'VariableNames', {'Depth','Sv','Pp','Sh_min','Sh_max','Differential','Brittleness'});
        resultsTable.Data = results;

        updateStressPlot();
        updateAdvancedPlots();
    end

    function updateStressPlot(~,~)
        cd = getappdata(fig,'ComputedData');
        if isempty(cd), return; end
        d = cd.depth;

        % Total stresses
        cla(axStress); hold(axStress,'on');
        legends = {};
        if chkVertical.Value
            plot(axStress, cd.S_v_MPa, d,'k-','LineWidth',2);
            legends{end+1} = 'Vertical Stress';
        end
        if chkShMin.Value
            plot(axStress, cd.S_h_min_MPa, d,'b--','LineWidth',2);
            legends{end+1} = 'Min Horizontal Stress';
        end
        if chkShMax.Value
            plot(axStress, cd.S_h_max_MPa, d,'r--','LineWidth',2);
            legends{end+1} = 'Max Horizontal Stress';
        end
        if chkPore.Value
            plot(axStress, cd.Pp_MPa, d,'g:','LineWidth',2);
            legends{end+1} = 'Pore Pressure';
        end
        hold(axStress,'off');
        set(axStress,'YDir','reverse');
        xlabel(axStress,'Stress (MPa)');
        ylabel(axStress,'Depth (m)');
        if ~isempty(legends), legend(axStress,legends,'Location','best'); end
        grid(axStress,'on');

        % Effective stresses
        cla(axEff);
        if chkEff.Value
            hold(axEff,'on');
            plot(axEff, cd.eff_S_v_MPa, d,'k-','LineWidth',2);
            plot(axEff, cd.eff_S_h_min_MPa, d,'b--','LineWidth',2);
            plot(axEff, cd.eff_S_h_max_MPa, d,'r--','LineWidth',2);
            hold(axEff,'off');
            legend(axEff,{'Eff. Vertical Stress','Eff. Min Horizontal Stress','Eff. Max Horizontal Stress'},'Location','best');
        end
        set(axEff,'YDir','reverse');
        xlabel(axEff,'Effective Stress (MPa)');
        ylabel(axEff,'Depth (m)');
        grid(axEff,'on');

        % Brittleness
        cla(axBrittle);
        scatter(axBrittle, cd.brittleness, d,36,cd.ShaleVolume,'filled');
        colormap(axBrittle, jet);
        cb = colorbar(axBrittle); cb.Label.String='ShaleVolume';
        set(axBrittle,'YDir','reverse');
        xlabel(axBrittle,'Brittleness (%)');
        ylabel(axBrittle,'Depth (m)');
        grid(axBrittle,'on');
    end

    function updateAdvancedPlots(~,~)
        cd = getappdata(fig,'ComputedData');
        if isempty(cd), return; end

        % 1) lambda vs mu
        cla(axLame);
        scatter(axLame, cd.lambda/1e9, cd.mu/1e9,36,cd.ShaleVolume,'filled');
        colormap(axLame,jet);
        cb1 = colorbar(axLame); cb1.Label.String='ShaleVolume';
        xlabel(axLame,'\lambda (GPa)'); ylabel(axLame,'\mu (GPa)');
        grid(axLame,'on');

        % 2) UCS vs Effective Stress
        UCS = 0.7 * cd.eff_S_v_MPa;
        horiz_eff = (cd.eff_S_h_min_MPa + cd.eff_S_h_max_MPa)/2;
        total_eff = (cd.eff_S_v_MPa + horiz_eff)/2;
        cla(axScatterUCS);
        scatter(axScatterUCS, UCS, total_eff,36,cd.ShaleVolume,'filled');
        colormap(axScatterUCS,jet);
        cb2 = colorbar(axScatterUCS); cb2.Label.String='ShaleVolume';
        xlabel(axScatterUCS,'UCS (MPa)');
        ylabel(axScatterUCS,'Effective Stress (MPa)');
        grid(axScatterUCS,'on');

        % 3) Horizontal Effective Stress vs Strain
        h_eff_strain = (horiz_eff*1e6 ./ cd.E_app)*100;
        cla(axScatterHoriz);
        scatter(axScatterHoriz, h_eff_strain, horiz_eff,36,cd.ShaleVolume,'filled');
        colormap(axScatterHoriz,jet);
        cb3 = colorbar(axScatterHoriz); cb3.Label.String='ShaleVolume';
        xlabel(axScatterHoriz,'Strain (%)');
        ylabel(axScatterHoriz,'Horizontal Effective Stress (MPa)');
        grid(axScatterHoriz,'on');

        % 4) Tensile Strength vs Porosity
        TensileStrength = UCS / 10;
        cla(axScatterVert);
        scatter(axScatterVert, TensileStrength, cd.Porosity,36,cd.ShaleVolume,'filled');
        colormap(axScatterVert,jet);
        cb4 = colorbar(axScatterVert); cb4.Label.String='ShaleVolume';
        xlabel(axScatterVert,'Tensile Strength (MPa)');
        ylabel(axScatterVert,'Porosity');
        grid(axScatterVert,'on');
    end

    function updateStatsTable(~,~)
        cd = getappdata(fig,'ComputedData');
        if isempty(cd)
            uialert(fig,'No computed data available.','Statistics Error');
            return;
        end
        fields = {'S_v_MPa','Pp_MPa','S_h_min_MPa','S_h_max_MPa','diff_MPa', ...
                  'eff_S_v_MPa','eff_S_h_min_MPa','eff_S_h_max_MPa','eff_diff_MPa','brittleness'};
        names  = {'Vert Stress','Pore Press','Min Hor Stress','Max Hor Stress','Diff Stress', ...
                  'Eff Vert Stress','Eff Min Hor','Eff Max Hor','Eff Diff','Brittleness'};
        statNames = {'Mean','Std Dev','Median','Min','Max','Std Error'};
        nParams = numel(fields);
        nStats  = numel(statNames);
        statsData = zeros(nStats,nParams);
        N = numel(cd.depth);
        for j = 1:nParams
            v = cd.(fields{j});
            statsData(1,j) = mean(v);
            statsData(2,j) = std(v);
            statsData(3,j) = median(v);
            statsData(4,j) = min(v);
            statsData(5,j) = max(v);
            statsData(6,j) = std(v)/sqrt(N);
        end
        set(statsTable,'Data',statsData, 'RowName',statNames,'ColumnName',names);

        x = 1:nParams;
        cla(axErrorStats);
        errorbar(axErrorStats, x, statsData(1,:), statsData(2,:), 'o-','LineWidth',2,'MarkerSize',8);
        set(axErrorStats,'XTick',x,'XTickLabel',names,'XTickLabelRotation',45);
        xlabel(axErrorStats,'Parameters');
        ylabel(axErrorStats,'Value');
        grid(axErrorStats,'on');
    end

    function updateSensitivityAnalysis(~,~)
        cd = getappdata(fig,'ComputedData');
        if isempty(cd)
            uialert(fig,'No computed data available.','Sensitivity Error');
            return;
        end
        fields = {'S_v_MPa','Pp_MPa','S_h_min_MPa','S_h_max_MPa','diff_MPa', ...
                  'eff_S_v_MPa','eff_S_h_min_MPa','eff_S_h_max_MPa','eff_diff_MPa','brittleness'};
        names  = {'Vert Stress','Pore Press','Min Hor Stress','Max Hor Stress','Diff Stress', ...
                  'Eff Vert Stress','Eff Min Hor','Eff Max Hor','Eff Diff','Brittleness'};
        nParams = numel(fields);
        sens = zeros(1,nParams);
        for j=1:nParams
            v = cd.(fields{j});
            m = mean(v); sd = std(v);
            sens(j) = (m~=0)*100*(sd/m);
        end
        cla(axSensitivity);
        bar(axSensitivity,1:nParams,sens,'FaceColor',[0.2 0.6 0.5]);
        set(axSensitivity,'XTick',1:nParams,'XTickLabel',names,'XTickLabelRotation',45);
        xlabel(axSensitivity,'Parameters');
        ylabel(axSensitivity,'Sensitivity (% C.V.)');
        grid(axSensitivity,'on');

        tbl = [names', num2cell(sens')];
        set(sensitivityTable,'Data',tbl,'ColumnName',{'Parameter','Sensitivity (%)'});
    end
end
