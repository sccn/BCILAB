function com = pop_topohistplot( EEG, typeplot, components, topotitle, rowcols, varargin);

com = '';
if nargin < 1
    help pop_topoplot;
    return;
end;
if nargin < 2
    typeplot = 1;
end;
if typeplot == 0 & isempty(EEG.icasphere)
    disp('Error: no ICA data for this set, first run ICA or Load AMICA components'); return;
end;
if isempty(EEG.chanlocs)
    disp('Error: cannot plot topography without channel location file'); return;
end;
if ~isfield(EEG.etc,'amica') || isempty(EEG.etc.amica)
    error('No AMICA solution found. You should first load AMICA components');return;
end




if nargin < 3
    % which set to save
    % -----------------
    if typeplot
        txtwhat2plot1 = 'Plotting ERP scalp maps at these latencies';
        txtwhat2plot2 = sprintf('(range: %d to %d ms, NaN -> empty):', ...
            round(EEG.xmin*1000), round(EEG.xmax*1000));
        editwhat2plot = [''];
    else
        txtwhat2plot1 = 'Component numbers';
        txtwhat2plot2 = '(negate index to invert component polarity; NaN -> empty subplot; Ex: -1 NaN 3)';
        editwhat2plot = ['1:' int2str(size(EEG.icaweights,1))];
    end;
    allmodels = 1:EEG.etc.amica.num_models;
    
    if ~isfield(EEG.etc.amica,'modnames')
        EEG.etc.amica.modnames{length(allmodels)} = '';
    end
    
    str = [num2str(allmodels(1)) fastif(~isempty(EEG.etc.amica.modnames{1}) ,[' - ' EEG.etc.amica.modnames{1}],'')];
    for i = 2:EEG.etc.amica.num_models
        str = [str '|' num2str(allmodels(i)) fastif(~isempty(EEG.etc.amica.modnames{i}) ,[' - ' EEG.etc.amica.modnames{i}],'')];
    end
    editwhat2select = str;
    
    if EEG.nbchan > 64,
        elecdef = ['''electrodes'', ''off'''];
    else,
        elecdef = ['''electrodes'', ''on'''];
    end;
    cb_enable = ['set(findobj(''parent'', gcbf, ''tag'', ''block'')    , ''enable'', ''on'');'];
    cb_disable = ['set(findobj(''parent'', gcbf, ''tag'', ''block'')    , ''enable'', ''off'');'];
    cb_block = ['if get(gcbo,''value''),' cb_enable ...
        'else,' cb_disable ...
        'end;'];
    cb_chosen = ['m = get(gcbo,''value''); set(findobj(''parent'',gcbf,''tag'',''title''),''string'',[EEG.setname '', Model '' num2str(m)])'];
    
    %EEG.etc.amica
    uilist = { { 'style'   'text'     'string'    txtwhat2plot1 } ...
        { 'style'   'edit'     'string'    editwhat2plot } ...
        { 'style'   'text'     'string'    txtwhat2plot2 } ...
        { } ...
        {'style' 'text' 'string' 'Select AMICA model'} ...
        {'style' 'popupmenu' 'string' editwhat2select 'value' 1 'callback' cb_chosen} ...
        { 'style'   'text'     'string'    'Plot title' } ...
        { 'style'   'edit'     'tag' 'title' 'string'    fastif(~isempty(EEG.setname), [EEG.setname ', Model 1'], '') } ...
        { 'style'   'text'     'string'    'Plot geometry (rows,col.); [] -> near square' } ...
        { 'style'   'edit'     'string'    '[]' } ...
        { 'style'   'text'     'string'    'Plot associated dipole(s) (if present)' } ...
        { 'style'   'checkbox' 'string'    '' } { } ...
        { } ...
        { 'style'   'text'     'string'    [ '-> Additional topoplot()' fastif(typeplot,'',' (and dipole)') ...
        ' options (see Help)' ] } ...
        { 'style'   'edit'     'string'    elecdef } ...
        {'style' 'checkbox' 'string' 'Show IC histograms' 'value' 1 'callback' cb_block} ...
        {'style' 'checkbox' 'tag' 'block' 'string' 'Use block analysis when determining histogram range' 'enable' 'on'}};
    uigeom = { [1.5 1] [1] [1] [1 1] [1.5 1] [1.5 1] [1.55 0.2 0.8] [1] [1] [1] [1] [1]};
    if typeplot
        uilist(9:11) = [];
        uigeom(6) = [];
    end;
    guititle = fastif( typeplot, 'Plot ERP scalp maps in 2-D -- pop_topoplot()', ...
        'Plot AMICA component scalp maps in 2-D -- pop_topohistplot()');
    
    result = inputgui( uigeom, uilist, 'pophelp(''pop_topoplot'')', guititle, [], 'normal');
    if isempty(result)
        return;
    end
    
    % reading first param
    % -------------------
    components   	     = eval( [ '[' result{1} ']' ] );
    if length(components) > EEG.nbchan
        tmpbut = questdlg2(...
            ['This involves drawing ' int2str(length(components)) ' plots. Continue ?'], ...
            '', 'Cancel', 'Yes', 'Yes');
        if strcmp(tmpbut, 'Cancel'), return; end;
    end;
    if isempty(components), error('Nothing to plot; enter parameter in first edit box'); end;
    
    % reading other params
    % --------------------
    h = result{2}; % chosen AMICA model
    topotitle   = result{3};
    rowcols     = eval( [ '[' result{4} ']' ] );
    options = [];
    if typeplot
        plotdip = 0;
        try, options      = eval( [ '{ ' result{5} ' }' ]);
        catch, error('Invalid scalp map options'); end;
    else
        plotdip     = result{5};
        try, options      = eval( [ '{ ' result{6} ' }' ]);
        catch, error('Invalid scalp map options'); end;
    end;
    if length(components) == 1,
        figure('paperpositionmode', 'auto'); curfig=gcf;
        try, icadefs;
            set(curfig, 'color', BACKCOLOR);
        catch, end;
    end;
    showhist = result{7};
    use_block = result{8};
    options = {options{:} 'showhist' showhist 'use_block' use_block 'model' h};
else
    if ~isempty(varargin) & isnumeric(varargin{1})
        plotdip = varargin{1};
        varargin = varargin(2:end);
    else
        plotdip = 0;
    end;
    options = varargin;
end;

% additional options
% ------------------
outoptions = { options{:} }; % for command
options    = { options{:} 'masksurf' 'on' };

% find maplimits
% --------------
maplimits = [];
for i=1:2:length(options)
    if isstr(options{i})
        if strcmpi(options{i}, 'maplimits')
            maplimits = options{i+1};
            options(i:i+1) = [];
            break;
        end;
    end;
end;

%--------------Added for IC histogram plotting--------------------------
for i=1:2:length(options)
    if isstr(options{i})
        if strcmpi(options{i}, 'showhist')
            showhist = options{i+1};
            options(i:i+1) = [];
            break;
        end;
    end;
end;

for i=1:2:length(options)
    if isstr(options{i})
        if strcmpi(options{i}, 'use_block')
            use_block = options{i+1};
            options(i:i+1) = [];
            break;
        end;
    end;
end;

for i=1:2:length(options)
    if isstr(options{i})
        if strcmpi(options{i}, 'model')
            h = options{i+1};
            options(i:i+1) = [];
            break;
        end;
    end;
end;
if ~exist('h','var')
    h = 1;
end
% ----------------------------------------------------------------------


nbgraph = size(components(:),1);
if ~exist('topotitle')
    topotitle = '';
end;
if ~exist('rowcols') | isempty(rowcols) | rowcols == 0
    rowcols(2) = ceil(sqrt(nbgraph));
    rowcols(1) = ceil(nbgraph/rowcols(2));
end;

SIZEBOX = 150;

fprintf('Plotting...\n');
if isempty(EEG.chanlocs)
    fprintf('Error: set has no channel location file\n');
    return;
end;

% Check if pop_topoplot input 'colorbar' was called, and don't send it to topoplot
loc = strmatch('colorbar', options(1:2:end), 'exact');
loc = loc*2-1;
if ~isempty(loc)
    colorbar_switch = strcmp('on',options{ loc+1 });
    options(loc:loc+1) = [];
else
    colorbar_switch = 1;
end

% determine the scale for plot of different times (same scales)
% -------------------------------------------------------------
if typeplot
    SIGTMP = reshape(EEG.data, EEG.nbchan, EEG.pnts, EEG.trials);
    pos = round( (components/1000-EEG.xmin)/(EEG.xmax-EEG.xmin) * (EEG.pnts-1))+1;
    nanpos = find(isnan(pos));
    pos(nanpos) = 1;
    SIGTMPAVG = mean(SIGTMP(:,pos,:),3);
    SIGTMPAVG(:, nanpos) = NaN;
    if isempty(maplimits)
        maxlim = max(SIGTMPAVG(:));
        minlim = min(SIGTMPAVG(:));
        maplimits = [ -max(maxlim, -minlim) max(maxlim, -minlim)];
    end;
else
    if isempty(maplimits)
        maplimits = 'absmax';
    end;
end;

if plotdip
    if strcmpi(EEG.dipfit.coordformat, 'CTF')
        disp('Cannot plot dipole on scalp map for CTF MEG data');
    end;
end;





%Calculate the winv matrix for the model chosen
A = pinv(EEG.etc.amica.W(:,:,h)*EEG.etc.amica.S);
%------------------------------------------------




% plot the graphs
% ---------------
counter = 1;
countobj = 1;
allobj = zeros(1,1000);
curfig = get(0, 'currentfigure');
if isfield(EEG, 'chaninfo'), options = { options{:} 'chaninfo' EEG.chaninfo }; end



if showhist
    
    alpha = EEG.etc.amica.alpha;
    sbeta = EEG.etc.amica.sbeta;
    rho = EEG.etc.amica.rho;
    mu = EEG.etc.amica.mu;
    nbgraph = length(components);
    mn = EEG.etc.amica.data_mean;
    x = EEG.data;
    
    models = EEG.etc.amica;
    if EEG.trials>1
        x = reshape(x,EEG.nbchan,EEG.trials*EEG.pnts);
        models.LLt = reshape(models.LLt,EEG.etc.amica.num_models,size(models.LLt,2)*size(models.LLt,3));
    end
    if max(size(models.LLt)) > 1
        indx = find(models.LLt(1,:)~=0);
        lim = find(indx>size(x,2));
        if ~isempty(lim)
            indx = indx(1:lim(1)-1);
        end
        x = x(:,indx);
        models.LLt = models.LLt(:,indx);
    end
    [nx,N] = size(x);
    n = size(models.c,1);
    B = 1000;
    num_blocks = floor(N/B);
    nbins = min(100,sqrt(N));
    [mval,mind] = max(models.LLt);
    n = size(models.c,1);
    s = zeros(n,N);
    disp('Sphering the data...')
    if use_block
        
        % sphere the data
        
        for k = 1:num_blocks
            xstrt = (k-1)*B + 1;
            if k < num_blocks
                xstp = k*B;
            else
                xstp = N;
            end
            for i = 1:nx
                x(i,xstrt:xstp) = x(i,xstrt:xstp) - mn(i);
                %x(i,xstrt:xstp) = s(1,xstrt:xstp) - mn(i);
            end
            x(1:n,xstrt:xstp) = models.S(1:n,:) * x(:,xstrt:xstp);
        end
        
        
        
        
        
        for k = 1:num_blocks
            xstrt = (k-1)*B + 1;
            if k < num_blocks
                xstp = k*B;
            else
                xstp = N;
            end
            s(components,xstrt:xstp) = models.W(components,:,h) * x(1:n,xstrt:xstp);
            for i = components
                s(i,xstrt:xstp) = s(i,xstrt:xstp) - models.c(i,h);
            end
        end
        
        
        
    else
        for i = 1:nx
            x(i,:) = x(i,:) - mn(i);
            %x(i,xstrt:xstp) = s(1,xstrt:xstp) - mn(i);
        end
        
        x(1:n,:) = models.S(1:n,:) * x(:,:);
        
        s(components,:) = models.W(components,:,h) * x(1:n,:);
        for i = components
            s(i,:) = s(i,:) - models.c(i,h);
        end
        
    end
    
end




for index = 1:size(components(:),1)
    if nbgraph > 1
        if mod(index, rowcols(1)*rowcols(2)) == 1
            
            if index> 1, figure(curfig); a = textsc(0.5, 0.05, topotitle); set(a, 'fontweight', 'bold'); end;
            curfig = figure('paperpositionmode', 'auto');
            pos = get(curfig,'Position');
            
            posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
            
            posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
            if showhist
                
                set(curfig,'Position', [posx posy  2*SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
            else
                set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
            end
            try, icadefs; set(curfig, 'color', BACKCOLOR); catch, end;
        end;
        
        if showhist
            curax = subplot( rowcols(1), 2*rowcols(2), 2*mod(index-1, rowcols(1)*rowcols(2))+1,'parent',curfig);
        else
            curax = subplot( rowcols(1), rowcols(2), mod(index-1, rowcols(1)*rowcols(2))+1,'parent',curfig);
        end
        set(curax, 'visible', 'off')
        if index == 1
            constantPosition = get(curax,'position');
        end
        
    else
        if showhist
            curax = subplot(1,2,1);
            set(curax,'visible','off')
        end
        
    end
    
    % add dipole location if present
    % ------------------------------
    dipoleplotted = 0;
    if plotdip && typeplot == 0
        if isfield(EEG, 'dipfit') & isfield(EEG.dipfit, 'model')
            if length(EEG.dipfit.model) >= index & ~strcmpi(EEG.dipfit.coordformat, 'CTF')
                %curpos = EEG.dipfit.model(components(index)).posxyz/EEG.dipfit.vol.r(end);
                curpos = EEG.dipfit.model(components(index)).posxyz;
                curmom = EEG.dipfit.model(components(index)).momxyz;
                try,
                    select = EEG.dipfit.model(components(index)).select;
                catch ,
                    select = 0;
                end;
                if ~isempty(curpos)
                    if strcmpi(EEG.dipfit.coordformat, 'MNI') % from MNI to sperical coordinates
                        transform = pinv( sph2spm );
                        tmpres = transform * [ curpos(1,:) 1 ]'; curpos(1,:) = tmpres(1:3);
                        tmpres = transform * [ curmom(1,:) 1 ]'; curmom(1,:) = tmpres(1:3);
                        try, tmpres = transform * [ curpos(2,:) 1 ]'; curpos(2,:) = tmpres(1:3); catch, end;
                        try, tmpres = transform * [ curmom(2,:) 1 ]'; curmom(2,:) = tmpres(1:3); catch, end;
                    end;
                    curpos = curpos / 85;
                    if size(curpos,1) > 1 && length(select) == 2
                        
                        dipole_index = find(strcmpi('dipole',options),1);
                        if  ~isempty(dipole_index) % if 'dipoles' is already defined in options{:}
                            options{dipole_index+1} = [ curpos(:,1:2) curmom(:,1:3) ];
                        else
                            options = { options{:} 'dipole' [ curpos(:,1:2) curmom(:,1:3) ] };
                        end
                        dipoleplotted = 1;
                    else
                        if any(curpos(1,:) ~= 0)
                            dipole_index = find(strcmpi('dipole',options),1);
                            if  ~isempty(dipole_index) % if 'dipoles' is already defined in options{:}
                                options{dipole_index+1} = [ curpos(1,1:2) curmom(1,1:3) ];
                            else
                                options = { options{:} 'dipole' [ curpos(1,1:2) curmom(1,1:3) ] };
                            end
                            dipoleplotted = 1;
                        end
                    end
                    
                end
                if nbgraph ~= 1
                    dipscale_index = find(strcmpi('dipscale',options),1);
                    if ~isempty(dipscale_index) % if 'dipscale' is already defined in options{:}
                        options{dipscale_index+1} = 0.6;
                    else
                        options = {  options{:} 'dipscale' 0.6 };
                    end
                end
                %options = { options{:} 'dipsphere' max(EEG.dipfit.vol.r) };
            end
        end
    end
    %curax = subplot( rowcols(1), 2*rowcols(2), 2*mod(index-1, rowcols(1)*rowcols(2))+1);
    % plot scalp map
    % --------------
    if index == 1
        addopt = { 'verbose', 'on' };
    else
        addopt = { 'verbose', 'off' };
    end;
    %fprintf('Printing to figure %d.\n',curfig);
    options = {  'maplimits' maplimits options{:} addopt{:} };
    if ~isnan(components(index))
        if typeplot
            if nbgraph > 1, axes(curax); end;
            tmpobj = topoplot( SIGTMPAVG(:,index), EEG.chanlocs, options{:});
            if nbgraph == 1,
                figure(curfig); if nbgraph > 1, axes(curax); end;
                title( [ 'Latency ' int2str(components(index)) ' ms from ' topotitle]);
            else
                figure(curfig); if nbgraph > 1, axes(curax); end;
                title([int2str(components(index)) ' ms'] );
            end;
        else
            if components(index) < 0
                figure(curfig);  if nbgraph > 1, axes(curax); end;
                tmpobj = topoplot( -A(:, -components(index)), EEG.chanlocs, options{:} );
            else
                figure(curfig);  if nbgraph > 1, axes(curax); end;
                tmpobj = topoplot( A(:, components(index)), EEG.chanlocs, options{:} );
            end;
            if nbgraph == 1, texttitle = [ 'IC ' int2str(components(index)) ' from ' topotitle];
            else             texttitle = ['' int2str(components(index))];
            end;
            if dipoleplotted, texttitle = [ texttitle ' (' num2str(EEG.dipfit.model(components(index)).rv*100,2) '%)']; end;
            figure(curfig);  if nbgraph > 1, axes(curax); end; title(texttitle);
        end;
        allobj(countobj:countobj+length(tmpobj)-1) = tmpobj;
        countobj = countobj+length(tmpobj);
        pos = get(curax,'position');
        
        set(curax,'position',[pos(1) pos(2) constantPosition(3) constantPosition(4)]);
        drawnow;
        axis square;
    else
        axis off
    end;
    
    % plot IC histogram
    % --------------------------
    
    if showhist
        
        
        
        
        if EEG.etc.amica.num_models > 1
            inds = find(mind == h);
            if ~isempty(inds)
                smn = mean(s(components(index),inds));
                [~,b] = hist(s(components(index),inds)-smn,nbins);
                
            else
                b = -10:0.01:10;
                
            end
        else
            [~,b] = hist(s(components(index),:),nbins);
        end
        
        
        
        if nbgraph > 1
            subplot(rowcols(1),2*rowcols(2),2*index,'parent',curfig);
        else
            subplot(1,2,2,'parent',curfig);
        end
        hold on
        sum = 0;
        for j = 1:size(alpha,1)
            mixtpdf = alpha(j,components(index),h)*sbeta(j,components(index),h)*(1/(2*gamma(1+1/rho(j,components(index),h))))*exp(-abs(sbeta(j,components(index),h)*(b-mu(j,components(index),h))).^rho(j,components(index),h));
            if ~isnan(mixtpdf)
                sum = mixtpdf + sum;
                plot(b,mixtpdf,'g');
            end
        end
        
        plot(b,sum,'b');
        
        drawnow;
    end
    % ---------------------------------------------
end

% Draw colorbar
if colorbar_switch
    if nbgraph == 1
        if ~isstr(maplimits)
            ColorbarHandle = cbar(0,0,[maplimits(1) maplimits(2)]);
        else
            ColorbarHandle = cbar(0,0,get(gca, 'clim'));
        end;
        pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
        if showhist %if IC histogram is also plotted
            set(ColorbarHandle,'position',[pos(1)+.02 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
        else
            set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.7 pos(4)-0.26]);
        end
    elseif ~isstr(maplimits)
        cbar('vert',0,[maplimits(1) maplimits(2)]);
    else cbar('vert',0,get(gca, 'clim'));
    end
    if ~typeplot    % Draw '+' and '-' instead of numbers for colorbar tick labels
        tmp = get(gca, 'ytick');
        set(gca, 'ytickmode', 'manual', 'yticklabelmode', 'manual', 'ytick', [tmp(1) tmp(end)], 'yticklabel', { '-' '+' });
    end
end

if nbgraph> 1,
    figure(curfig); a = textsc(0.5, 0.05, topotitle);
    set(a, 'fontweight', 'bold');
end;
if nbgraph== 1,
    com = 'figure;';
end;
set(allobj(1:countobj-1), 'visible', 'on');

figure(curfig);
axcopy(curfig, 'set(gcf, ''''units'''', ''''pixels''''); postmp = get(gcf, ''''position''''); set(gcf, ''''position'''', [postmp(1) postmp(2) 560 420]); clear postmp;');

com = [com sprintf('pop_topohistplot(%s,%d, %s);', ...
    inputname(1), typeplot, vararg2str({components topotitle rowcols plotdip outoptions{:} }))];
disp('Done.');
return;
