function [mtable, buttons] = createTable(varargin)
% createTable - create nice-looking table based on javax.swing.JTable
%
% Syntax:
%    [mtable, buttons] = createTable (pnContainer, headers, data, buttonsFlag, 'PropName',PropValue, ...)
%
% Input Parameters:
%    pnContainer - optional handle to container uipanel or figure. If empty/unsupplied then current figure will be used
%    headers     - optional cell array of column header strings. If unsupplied then = {'A','B','C'}
%    data        - optional vector/matrix (either scalar or cell array) of data values
%    buttonsFlag - optional flag indicating whether to create the table-manipulation buttons. Default = true
%    'PropName',PropValue - 
%                  optional list of property pairs (e.g., 'AutoResizeMode',4,'Editable',false,'Position',[.1,.1,.5,.5])
%                  Note: PropName is either an mtable property ('Visible','Editable','Position','DataChangedCallback',...)
%                        or a Javax.swing.JTable property ('ShowGrid','Name','RowHeight',...)
%                        or a javax.swing.table.JTableHeader property ('ResizingAllowed',ReorderingAllowed',...).
%                        Abbreviated PropNames are unsupported for mtable properties (which are few) - only for JTable
%                  Note: All optional parameters of createTable may be specified using PropName/PropValue pairs,
%                        case-insensitive, in whichever order (see the bottom example below):
%                        - 'Container'                            (default=gcf)
%                        - 'Headers' (default={'A','B',...} based on data size)
%                        - 'Data'                                 (default=[])
%                        - 'Buttons' => 'on','off',1,0,true,false (default='on')
%                        - 'ImageColumns' => index or name        (default={})
%                        - 'ImageTooltipHeight'                   (default=300)
%                        - 'uicontextmenu'                        (default=[])
%
% Output parameters:
%    mtable      - handle to mtable object (a Matlab object)
%    buttons     - handles to table manipulation buttons: [<appendRow> <insertRow> <deleteRow> <deleteAll> <printAll>]
%
% Examples:
%    [mtable, buttons] = createTable;
%    [mtable, buttons] = createTable(gcf,'column name');
%    mtable = createTable([],{'a','b','c','d'},{false,1.3,uint16(45),'ert'; true,pi,uint16(-4),'defrgt'})
%    mtable = createTable([],{'a','b','c','d'},magic(4),false,'AutoResizeMode',javax.swing.JTable.AUTO_RESIZE_ALL_COLUMNS)
%    mtable = createTable([],{'rads','sin','cos'},[pi,sin(pi),cos(pi)],'SelectionMode',javax.swing.ListSelectionModel.SINGLE_INTERVAL_SELECTION)
%    mtable = createTable('Data',magic(3), 'Headers',{'a','b','c'}, 'Buttons','off', 'Container',gcf);
%
% Usage:
%    The table automatically resizes to fill the pnContainer (you may modify this via the 'Position' property).
%    The table automatically sets the columns' cell editor and renderer based on the supplied data. Logical values are
%       given a checkbox, strings are left-aligned (numbers are right-aligned). You can always override the defaults.
%    You can change column widths by dragging the column borders on the header row.
%    You can sort columns by clicking the column header (once to sort descending, once again to sort ascending and once
%       more for the unsorted view). Sorting multiple columns is done by control-clicking all relevant columns (the
%       sorting icon is decreased in size for each additional minor sort col).
%    You can copy/paste any consecutive region of table cells, just as in Excel. You can select entire rows or columns
%       by right-clicking their header. You can also paste Excel data directly, with Ctrl-Shift-V (or use the context
%       menu by right-clicking) at the target table cell.
%    For additional tips about how to set multiple aspects of the table, refer to:
%       <a href="http://java.sun.com/docs/books/tutorial/uiswing/components/table.html">http://java.sun.com/docs/books/tutorial/uiswing/components/table.html</a>
%
% Programming tips/cues/examples:
%    mtable = createTable(...)
%    jtable = mtable.getTable;
%    mtable.setVisible(false);
%    mtable.setCheckBoxEditor(1);  % Set first column to a checkbox (see Note 2 below)
%    cb = javax.swing.JComboBox({'First','Last'}); cb.setEditable(true);  % prepare an editable drop-down CellEditor
%    editor = javax.swing.DefaultCellEditor(cb);
%    jtable.getColumnModel.getColumn(1).setCellEditor(editor);  % assign this editor to second column (see Note 2)
%    jtable.getColumnModel.getColumn(0).setMaxWidth(20);  % Limit width of first (checkbox) column (see Note 2)
%    mtable.setEditable(0,false);  % Disable editing first column (see note 2 below)
%    renderer = javax.swing.table.DefaultTableCellRenderer;  % or: renderer = jtable.getColumnModel.getColumn(1).getCellRenderer
%    renderer.setHorizontalAlignment(javax.swing.SwingConstants.RIGHT);  % useful for numbers rendered as strings e.g.: jtable.setValueAt(sprintf('%.1f',pi,rowIdx,colIdx))
%    jtable.getColumnModel.getColumn(1).setCellRenderer(renderer);  % right-align second column (see note 2)
%    data = cell(mtable.getData);  % a cell matrix (mtable.getData is a java.lang.Object[][] object, using base-1 indexing)
%    data = mtable.getTableModel.getDataVector;  % a java.util.Vector object ([[false, 1.3, 45, ert], [true, 3.14,...]])
%    jtable.setValueAt(value,rowIdx,colIdx);  % 0-based Idx - see Note 2 below
%    jtable.getModel.addRow({true, pi, int16(45), 'test'});  % appends a row to the bottom of the table
%    mtable.DataChangedCallback = [];  % used to temporarily disable data-change callbacks
%    mtable.DataChangedCallback = @myDataChange_Callback;  % myDataChange_Callback is a Matlab function
%
%    % Sample dataChange_Callback function
%    function dataChange_Callback(mtable, eventdata)
%       if ~ishandle(mtable),  return;  end
%          % Prevent re-entry here if the callback is not thread-safe - see Note 3 below
%       eventDetails = eventdata.getEvent;
%       modifiedColIdx = eventDetails.getColumn;
%       modifiedRowIdx = eventDetails.getFirstRow;
%       if modifiedColIdx>=0 && modifiedRowIdx>=0
%          data = mtable.getData;
%          newValue = data(modifiedRowIdx+1,modifiedColIdx+1);  % see Note 2 below
%          switch modifiedColIdx
%            case ...
%          end
%       end
%
% Notes:
%    1. Some (very few) JTable features are inconsistent or unavailable in different jave versions. Type 
%       '<a href="matlab:version -java">version -java</a>' at the command prompt to see your specific java version.
%    2. Note that java uses 0-based indexing, while Matlab is 1-based. The returned mtable parameter is a Matlab object
%       (so use 1-base), while mtable.getXXX returns java objects (0-based). jtable above is an example of a java object.
%    3. Modifying mtable.DataChangedCallback within the callback doesn't work - you need to use some global flag/mutex
%    4. The <Print> button uses Excel to parse and print the table
%    5. Due to Matlab limitations (specifically, of uitable/UitablePeer) the table is created as a direct child of
%       the container figure (although it is visually positioned within pnContainer)
%    6. To enable sorting functionality, the attached TableSorter.jar file must be located in the java classpath.
%       See the Matlab documentation for <a href="matlab:doc javaclasspath">javaclasspath</a>. Note that using 
%       javaaddpath(...) to set the path has a nasty side-effect (at least since Matlab 7.2) of clearing all globals!
%       An alternative is to place the pathname for TableSorter.jar in the <a href="matlab:which classpath.txt">classpath.txt</a> file
%
% Known issues/limitations:
%    - Column alignment not preserved during Print
%    - Print fails if Excel unavailable (maybe directly print tab-separated text data)
%    - Unable to add/delete rows or to print via context menu (right-click)
%    - Table is created as a direct child of figure, not pnContainer (see Note 5 above)
%
% Bugs and suggestions:
%    Please send to Yair Altman (altmany at gmail dot com)
%
% See also:
%    uitable, java, javaclasspath
%
% Release history:
%    1.0 2007-03-09: initial version
%    1.1 2007-03-22: fixed selected row# on deletion of bottom row, main comment, missing option; added header tooltip
%    1.2 2007-03-25: fixed row deletion when entire row selected
%    1.3 2007-05-07: fix for Matlab 7.0.4 as per Sebastian Hölz; fix default headers per Selig Sechzer; full P-V parameters processing
%    1.4 2009-05-25: support for image cells, column reordering/resizing, uicontextmenu
%    1.5 2013-06-27: fix for HG2
%    1.6 2015-01-07: fixes for HG2 and a couple of input params edge-cases

% License to use and modify this code is granted freely to all interested, as long as the original author is
% referenced and attributed as such. The original author maintains the right to be solely associated with this work.

% Programmed and Copyright by Yair M. Altman: altmany(at)gmail.com
% $Revision: 1.6 $  $Date: 2015/01/07 21:52:53 $

  %try
      % Ensure that java swing is enabled...
      if ~usejava('swing')
          error('createTable:NeedSwing','Java tables require Java Swing.');
      end

      % Process optional arguments
      paramsStruct = processArgs(varargin{:});

      hContainer = handle(paramsStruct.container);
      if isa(hContainer,'figure') || isa(hContainer,'matlab.ui.Figure')
          pnContainerPos = getpixelposition(paramsStruct.container,0);  % Fix for Matlab 7.0.4 as per Sebastian Hölz
          pnContainerPos(1:2) = 0;
      else
          pnContainerPos = getpixelposition(paramsStruct.container,1);  % Fix for Matlab 7.0.4 as per Sebastian Hölz
      end

      % Get handle to parent figure
      hFig = ancestor(paramsStruct.container,'figure');

      % Determine whether table manipulation buttons are requested
      if paramsStruct.buttons
          margins = [1,30,0,-30];  % With buttons
      else
          margins = [1,1,0,0];   % No buttons
      end

      % Get the uitable's required position within the container
      tablePosition = pnContainerPos + margins;    % Relative to the figure

      % Start with dummy data, just so that uitable can be initialized (or use supplied data, if available)
      if isempty(paramsStruct.data)
          numRows = 0;
          numCols = length(paramsStruct.headers);
          paramsStruct.data = zeros(1,numCols);
      else
          numRows = size(paramsStruct.data,1);
      end

      % Create a sortable uitable within the container
      try
          % use the old uitable (Matlab R2008+)
          mtable = uitable('v0', hFig, 'position',tablePosition, 'Data',paramsStruct.data, 'ColumnNames',paramsStruct.headers);
          warning off MATLAB:uitable:DeprecatedFunction  % TODO: remove deprecation warning
      catch
          mtable = uitable(hFig, 'position',tablePosition, 'Data',paramsStruct.data, 'ColumnNames',paramsStruct.headers);
      end
      mtable.setNumRows(numRows);
      set(mtable,'units','normalized');  % this will resize the table whenever its container is resized

      % jtable is the underlying java JTable - access to lots more functionality...
      % Note: actually, jtable is a com.mathworks.hg.peer.UitablePeer$PeerSpreadsheetTable object, but this extends
      % ^^^^  javax.swing.JTable, so for all practical purposes you may use it as a JTable
      jtable = mtable.getTable;

      % Fix for JTable focus bug : see http://bugs.sun.com/bugdatabase/view_bug.do;:WuuT?bug_id=4709394
      % Taken from: http://xtargets.com/snippets/posts/show/37
      jtable.putClientProperty('terminateEditOnFocusLost', java.lang.Boolean.TRUE);

      % We want to use sorter, not data model...
      % unfortunately, UitablePeer expects DefaultTableModel (not TableSorter) so we need a modified UitablePeer class
      % however, UitablePeer is a Matlab class, so instead let's use a modified TableSorter and attach it to the Model
      %sorter = com.mathworks.toolbox.dasstudio.util.TableSorter;  % Failed attempt...
      %sorter = com.mathworks.mwswing.DefaultSortableTable;        % ...another failed attempt...
      if ~isempty(which('TableSorter'))
          % Add TableSorter as TableModel listener
          sorter = TableSorter(jtable.getModel());  %(table.getTableModel);
          %tablePeer = UitablePeer(sorter);  % This is not accepted by UitablePeer... - see comment above
          jtable.setModel(sorter);
          sorter.setTableHeader(jtable.getTableHeader());

          % Set the header tooltip (with sorting instructions)
          jtable.getTableHeader.setToolTipText('<html>&nbsp;<b>Click</b> to sort up; <b>Shift-click</b> to sort down<br>&nbsp;<b>Ctrl-click</b> (or <b>Ctrl-Shift-click</b>) to sort secondary&nbsp;<br>&nbsp;<b>Click again</b> to change sort direction<br>&nbsp;<b>Click a third time</b> to return to unsorted view<br>&nbsp;<b>Right-click</b> to select entire column</html>');
      else
          % Set the header tooltip (no sorting instructions...)
          jtable.getTableHeader.setToolTipText('<html>&nbsp;<b>Click</b> to select entire column<br>&nbsp;<b>Ctrl-click</b> (or <b>Shift-click</b>) to select multiple columns&nbsp;</html>');
      end

      % Store the uitable's handle within the container's userdata, for later use
      set(paramsStruct.container,'userdata',[get(paramsStruct.container,'userdata'), mtable]);  % add to parent userdata, so we have a handle for deletion

      % Enable multiple row selection, auto-column resize, and auto-scrollbars
      scroll = mtable.TableScrollPane;
      scroll.setVerticalScrollBarPolicy(scroll.VERTICAL_SCROLLBAR_AS_NEEDED);
      scroll.setHorizontalScrollBarPolicy(scroll.HORIZONTAL_SCROLLBAR_AS_NEEDED);
      jtable.setSelectionMode(javax.swing.ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);
      jtable.setAutoResizeMode(jtable.AUTO_RESIZE_SUBSEQUENT_COLUMNS)

      % Set the jtable name based on the containing panel's tag
      basicTagName = get(paramsStruct.container,'tag');
      jtable.setName([basicTagName 'Table']);

      % Move the selection to first table cell (if any data available)
      if (jtable.getRowCount > 0)
          jtable.changeSelection(0,0,false,false);
      end

      % Process optional args
      processParams(paramsStruct,mtable,jtable);

      % Create table manipulation buttons
      if paramsStruct.buttons
          buttons = createManipulationButtons(paramsStruct.container,mtable);
      else
          buttons = [];
      end
  %catch
      % Insert your code here
      %handleError;
  %end

%% Process optional arguments
function paramsStruct = processArgs(varargin)

    % Get the properties in either direct or P-V format
    [regParams, pvPairs] = parseparams(varargin);

    % Fix args in case of P-V mismatch
    if mod(numel(pvPairs),2)
        regParams{end+1} = pvPairs{1};
        pvPairs(1) = [];
    end

    % Now process the optional P-V params
    try
        % Initialize
        paramName = [];
        paramsStruct = [];
        paramsStruct.container = [];
        paramsStruct.headers = {'A','B','C'};  % 3 columns by default
        paramsStruct.data = {};
        paramsStruct.buttons = true;
        paramsStruct.imagecolumns = {};
        paramsStruct.imagetooltipheight = 300;  % 300px by default (max)
        paramsStruct.extra = {};

        % Parse the regular (non-named) params in recption order
        if length(regParams)>0,  paramsStruct.container = regParams{1};  end  %#ok
        if length(regParams)>1,  paramsStruct.headers   = regParams{2};  end
        if length(regParams)>2,  paramsStruct.data      = regParams{3};  end
        if length(regParams)>3,  paramsStruct.buttons   = regParams{4};  end

        % Parse the optional param PV pairs
        supportedArgs = {'container','headers','data','buttons','imagecolumns','imagetooltipheight'};
        while ~isempty(pvPairs)

            % Ensure basic format is valid
            paramName = '';
            if ~ischar(pvPairs{1})
                error('YMA:createTable:invalidProperty','Invalid property passed to createTable');
            elseif length(pvPairs) == 1
                error('YMA:createTable:noPropertyValue',['No value specified for property ''' pvPairs{1} '''']);
            end

            % Process parameter values
            paramName  = pvPairs{1};
            paramValue = pvPairs{2};
            pvPairs(1:2) = [];
            if any(strncmpi(paramName,supportedArgs,length(paramName)))
                paramsStruct.(lower(paramName)) = paramValue;
            else
                paramsStruct.extra = {paramsStruct.extra{:} paramName paramValue};
            end
        end  % loop pvPairs

        % Create a panel spanning entire figure area, if container handle was not supplied
        if isempty(paramsStruct.container) || ~ishandle(paramsStruct.container)
            paramsStruct.container = uipanel('parent',gcf,'tag','TablePanel');
        end

        % Set default header names, if not supplied
        if isempty(paramsStruct.headers)
            if isempty(paramsStruct.data)
                paramsStruct.headers = {' '};
            else
                paramsStruct.headers = cellstr(char('A'-1+(1:size(paramsStruct.data,2))'))';
            end
        elseif ischar(paramsStruct.headers)
            paramsStruct.headers = {paramsStruct.headers};
        end

        % Convert data to cell-format (if not so already)
        if ~iscell(paramsStruct.data)
            numCols = size(paramsStruct.data,2);
            paramsStruct.data = mat2cell(paramsStruct.data,ones(1,size(paramsStruct.data,1)),ones(1,numCols));
        end

        % Ensure a logical-convertible buttons flag
        if ischar(paramsStruct.buttons)
            switch lower(paramsStruct.buttons)
                case 'on',  paramsStruct.buttons = true;
                case 'off', paramsStruct.buttons = false;
                otherwise
                    error('YMA:createTable:invalidProperty','Invalid buttons property value: must be ''on'', ''off'', 1, 0, true or false');
            end
        elseif isempty(paramsStruct.buttons) || ~(isnumeric(paramsStruct.buttons) || islogical(paramsStruct.buttons))
            error('YMA:createTable:invalidProperty','Invalid buttons property value: must be ''on'', ''off'', 1, 0, true or false');
        end
    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        error('YMA:createTable:invalidProperty',['Error setting createTable property' paramName ':' char(10) lasterr]);
    end
%end  % processArgs


%% Process optional arguments on the newly-created table object
function processParams(paramsStruct,mtable,jtable)
    try
        % Process regular extra parameters
        paramName = '';
        th = jtable.getTableHeader;
        container = get(mtable,'uicontainer');
        for argIdx = 1 : 2 : length(paramsStruct.extra)
            if argIdx<2
                % We need this pause to let java complete all table rendering
                % TODO: We should really use calls to awtinvoke() instead, though...
                pause(0.05);
            end
            if (length(paramsStruct.extra) > argIdx)   % ensure the arg value is there...
                paramsStruct.extra{argIdx}(1) = upper(paramsStruct.extra{argIdx}(1));  % property names always start with capital letters...
                paramName  = paramsStruct.extra{argIdx};
                paramValue = paramsStruct.extra{argIdx+1};
                propMethodName = ['set' paramName];
                
                % First try to modify the container
                try
                    set(container, paramName, paramValue);
                catch
                    try % if ismethod(mtable,propMethodName)
                        % No good, so try the mtable...
                        set(mtable, paramName, paramValue);
                    catch %elseif ismethod(jtable,propMethodName)
                        try
                            % sometimes set(t,x,y) failes but t.setX(y) is ok...
                            javaMethod(propMethodName, mtable, paramValue);
                        catch
                            try
                                % Try to modify the underlying JTable itself
                                if isprop(jtable, paramName)
                                    set(jtable, paramName, paramValue);
                                else
                                    error('noSuchProp');
                                end
                            catch
                                try
                                    javaMethod(propMethodName, jtable, paramValue);
                                catch
                                    try
                                        % Try to modify the table header...
                                        set(th, paramName, paramValue);
                                    catch
                                        javaMethod(propMethodName, th, paramValue);
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end  % for argIdx

        % Process requested image columns
        if ~isempty(which('ImageCellRenderer')) && ~isempty(paramsStruct.imagecolumns)
            if ischar(paramsStruct.imagecolumns)
                % Maybe a header name?
                jtable.getColumn(paramsStruct.imagecolumns).setCellRenderer(ImageCellRenderer(paramsStruct.imagetooltipheight));
            elseif iscellstr(paramsStruct.imagecolumns)
                % Cell array of header names
                for argIdx = 1 : length(paramsStruct.imagecolumns)
                    jtable.getColumn(paramsStruct.imagecolumns{argIdx}).setCellRenderer(ImageCellRenderer(paramsStruct.imagetooltipheight));
                    drawnow;
                end
            else
                % Try to treat as a numeric index array
                for argIdx = 1 : length(paramsStruct.imagecolumns)
                    colIdx = paramsStruct.imagecolumns(argIdx) - 1;  % assume 1-based indexing
                    %jtable.setEditable(colIdx,0);  % images are editable!!!
                    jtable.getColumnModel.getColumn(colIdx).setCellRenderer(ImageCellRenderer(paramsStruct.imagetooltipheight));
                    drawnow;
                end
            end
             drawnow;
        elseif ~isempty(paramsStruct.imagecolumns)  % i.e., missing Renderer
            warning('YMA:createTable:missingJavaClass','Cannot set image columns: ImageCellRenderer.class is missing from the Java class path');
        end
        jtable.repaint;

        % Process UIContextMenu
        try cm = get(container,'uicontextmenu'); catch, cm=[]; end  % fails in HG2
        if ~isempty(cm)
            popupMenu = jtable.getRowHeaderPopupMenu;
            %popupMenu.list;
            popupMenu.removeAll; drawnow; pause(0.1);
            cmChildren = get(cm,'child');
            itemNum = 0;
            for cmChildIdx = length(cmChildren) : -1 : 1
                if itemNum == 6
                    % add 2 hidden separators which will be removed by the Matlab mouse listener...
                    popupMenu.addSeparator;
                    popupMenu.addSeparator;
                    popupMenu.getComponent(5).setVisible(0);
                    popupMenu.getComponent(6).setVisible(0);
                    itemNum = 8;
                end
                % Add a possible separator
                if strcmpi(get(cmChildren(cmChildIdx),'Separator'),'on')
                    popupMenu.addSeparator;
                    itemNum = itemNum + 1;
                end
                if itemNum == 6
                    % add 2 hidden separators which will be removed by the Matlab mouse listener...
                    popupMenu.addSeparator;
                    popupMenu.addSeparator;
                    popupMenu.getComponent(5).setVisible(0);
                    popupMenu.getComponent(6).setVisible(0);
                    itemNum = 8;
                end
                % Add the main menu item
                jMenuItem = javax.swing.JMenuItem(get(cmChildren(cmChildIdx),'Label'));
                set(jMenuItem,'ActionPerformedCallback',get(cmChildren(cmChildIdx),'Callback'));
                popupMenu.add(jMenuItem);
                itemNum = itemNum + 1;
            end
            for extraIdx = itemNum+1 : 7
                popupMenu.addSeparator;
                popupMenu.getComponent(extraIdx-1).setVisible(0);
            end
            drawnow;
        end

    catch
        if ~isempty(paramName),  paramName = [' ''' paramName ''''];  end
        error('YMA:createTable:invalidProperty',['Error setting createTable property' paramName ':' char(10) lasterr]);
    end
%end  % processParams


%% --- Executes on button press in btInsert.
function buttons = createManipulationButtons(pnContainer, mtable)
% pnContainer  handle to container uipanel
% mtable       handle to mtable (Matlab) object
  %try
      btAppendRow = uicontrol('tag','btTableAppendRow', 'callback',@btTableAppendRow_Callback, 'position',  [10,5,60,20], 'string','Append',     'parent',pnContainer, 'userdata',mtable);
      btInsertRow = uicontrol('tag','btTableInsertRow', 'callback',@btTableInsertRow_Callback, 'position',  [80,5,60,20], 'string','Insert',     'parent',pnContainer, 'userdata',mtable);
      btDeleteRow = uicontrol('tag','btTableDeleteRow', 'callback',@btTableDeleteRow_Callback, 'position', [150,5,60,20], 'string','Delete',     'parent',pnContainer, 'userdata',mtable);
      btDeleteAll = uicontrol('tag','btTableDeleteAll', 'callback',@btTableDeleteAll_Callback, 'position', [220,5,60,20], 'string','Delete All', 'parent',pnContainer, 'userdata',mtable);
      btPrintAll  = uicontrol('tag','btTablePrintAll',  'callback',@btTablePrintAll_Callback,  'position', [290,5,60,20], 'string','Print',      'parent',pnContainer, 'userdata',mtable);
      buttons = [btInsertRow btAppendRow btDeleteRow btDeleteAll btPrintAll];
      if mtable.getNumRows < 1
          setVisibility(pnContainer, 'off');
      end
  %catch
      % Insert your code here
      %handleError;
  %end


%% --- Executes on button press in btTableInsert.
% Insert a new row immediately before the current row
function btTableInsertRow_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to btTableInsertRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  try
      mtable = get(hObject,'userdata');
      jtable = mtable.getTable;

      % Stop any current editing
      stopEditing(jtable);

      % Insert the new row immediately before the current row
      newRowData = cell(1,mtable.getNumColumns);  % empty data
      mtable.getTableModel.insertRow(max(0,jtable.getSelectedRow), newRowData);
      jtable.repaint;
  catch
      % Insert your code here
      handleError;
  end


%% --- Executes on button press in btTableAppend.
% Insert a new row as the last row in the table
function btTableAppendRow_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to btTableAppendRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  %try
      mtable = get(hObject,'userdata');
      jtable = mtable.getTable;

      % Stop any current editing
      stopEditing(jtable);

      % Add a new row at the bottom of the data table
      newRowData = cell(1,mtable.getNumColumns);  % empty data
      mtable.getTableModel.addRow(newRowData);
      jtable.repaint;

      % Move the selection to Column A of this new row
      jtable.changeSelection(jtable.getRowCount-1,0,false,false);

      % There must be at least one table row now, so display the table in any case
      try mtable.setVisible(true); catch, end
      setVisibility(hObject, 'on');
  %catch
      % Insert your code here
      %handleError;
  %end


%% --- Executes on button press in btTableDelete.
% If there are any rows displayed, then delete the currently-selected row
function btTableDeleteRow_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to btTableDeleteRow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  %try
      mtable = get(hObject,'userdata');
      jtable = mtable.getTable;

      % Stop any current editing
      stopEditing(jtable);

      % If there are any rows displayed, then delete the currently-selected row
      rowCount = jtable.getRowCount;
      if (rowCount > 0)  % might be==0 during slow processing & user double-click
          currentRow = max(0,jtable.getSelectedRow);
          currentCol = max(0,jtable.getSelectedColumn);
          mtable.getTableModel.removeRow(currentRow);
          if currentRow >= rowCount-1
              jtable.changeSelection(currentRow-1, currentCol, false, false);
          elseif jtable.getSelectedRow < 0
              jtable.changeSelection(currentRow, currentCol, false, false);
          end
      end
      if (jtable.getRowCount <= 0)
          %table.setVisible(false);
          setVisibility(hObject, 'off');
      end
      jtable.repaint;
  %catch
      % Insert your code here
      %handleError;
  %end


%% --- Executes on button press in btTableDeleteAll.
% Deletes all table rows
function btTableDeleteAll_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to btTableDeleteAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  %try
      mtable = get(hObject,'userdata');
      jtable = mtable.getTable;

      % Stop any current editing
      stopEditing(jtable);

      % Delete all table rows
      mtable.setNumRows(0);

      % Hide irrelevant controls
      %mtable.setVisible(false);
      setVisibility(hObject, 'off');
  %catch
      % Insert your code here
      %handleError;
  %end


%% --- Executes on button press in btTablePrint.
% Prints the table via Excel
function btTablePrintAll_Callback(hObject, eventdata, handles)  %#ok
% hObject    handle to btTablePrintAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
  persistent hExcel
  try
      mtable = get(hObject,'userdata');

      % Try to open an Excel COM server
      % Note: http://msdn.microsoft.com/library/default.asp?url=/library/en-us/odc_2003_ta/html/odc_landoffice03_vba.asp
      try
          % try to reuse an existing (pre-opened) COM server
          % If we can't access the ActiveX parent, it means it's closed
          parent = hExcel.parent;  %#ok
      catch
          try
              % Try to reuse an existing Excel COM server if possible
              hExcel = actxGetRunningServer('excel.application');
          catch
              hExcel = actxserver('excel.application');
          end
      end

      % Try to open the requested document
      hExcel.Workbooks.Add;

      % Format field headers
      headers = cell(mtable.getColumnNames)';
      if ~isempty(headers)
          hExcel.Range(['A1:' n2a(length(headers)) '1']).Select;
          hExcel.Selection.Value = headers;
          hExcel.Selection.Font.Bold = true;                   % bold
          hExcel.Selection.Font.Color = hex2dec('0000FF');     % Red
          hExcel.Selection.Border.Item(4).Weight = 3;          % underline
          hExcel.Selection.Cells.HorizontalAlignment = -4108;  % =xlCenter
      end

      % Set the data from the table
      data = cell(mtable.data);
      if ~isempty(headers)
          hExcel.Range(['A2:' n2a(size(data,2)) num2str(1+size(data,1))]).Select;
          hExcel.Selection.Value = data;
      end

      % TODO: Change checkbox fields to boolean (TRUE/empty)

      % Other formats
      %hExcel.Cells.HorizontalAlignment = -4108;  % =xlCenter  % TODO: preserve original jtable column alignment
      hExcel.Cells.EntireColumn.AutoFit;
      hExcel.ActiveSheet.DisplayRightToLeft = false;
      set(hExcel.ActiveSheet.PageSetup, 'LeftMargin',  hExcel.InchesToPoints(0.1), ...
                                        'RightMargin', hExcel.InchesToPoints(0.1), ...
                                        'HeaderMargin',hExcel.InchesToPoints(0), ...
                                        'FooterMargin',hExcel.InchesToPoints(0.1), ...
                                        'TopMargin',120, ...
                                        'BottomMargin',36, ...
                                        'FitToPagesWide',1, ...
                                        'FitToPagesTall',1, ...
                                        'Orientation','xlPortrait', ...
                                        'LeftHeader','&D  &T', ...
                                        'CenterHeader',char(mtable.getTable.getName), ...
                                        'RightHeader','&G');

      % Send to printer
      hExcel.ActiveWindow.SelectedSheets.PrintOut;

      % Close the workbook
      invoke(hExcel.ActiveWindow,'close',false);
  catch
      % Insert your code here
      %handleError;
      err = lasterror;
      try  invoke(hExcel.ActiveWindow,'close',false);  catch  end;  % just in case of a printing error
      rethrow(err);
  end


%% --- Convert col # format to 'A','B','C','AA',... format
% Thanks Brett Shoelson, via CSSM
function colStr = n2a(c)
  t = [floor((c-1)/26)+64, rem(c-1,26)+65];
  if (t(1)<65), t(1) = []; end
  colStr = char(t);


%% --- Executes on button press in btInsert.
function stopEditing(jtable)
  %try
      component = jtable.getEditorComponent;
      if ~isempty(component)
          event = javax.swing.event.ChangeEvent(component);
          jtable.editingStopped(event);
      end
  %catch
      % Insert your code here
      %handleError;
  %end

  
%% --- Utility function to set visibility of row manipulation buttons
function setVisibility(hObject, enableStr)
% hObject    handle to some element within the figure
% enableStr  'on' or 'off'
  %try
      hParent = ancestor(hObject,'figure');
      set(findall(hParent,'tag','btTableInsertRow'),'enable',enableStr);
      set(findall(hParent,'tag','btTableDeleteRow'),'enable',enableStr);
      set(findall(hParent,'tag','btTableDeleteAll'),'enable',enableStr);
      set(findall(hParent,'tag','btTablePrintAll'), 'enable',enableStr);
  %catch
      % Insert your code here
      %handleError;
  %end
