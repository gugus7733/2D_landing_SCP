classdef FigureManager < handle
    %% FIGUREMANAGER – reusable batch-plot helper
    properties
        cfg struct
        handles struct = struct()
        idx double = 0
        t0stamp char
    end

    %% PUBLIC =============================================================
    methods
        function obj = FigureManager(varargin)
            obj.cfg = obj.defaultCfg();
            obj.cfg = obj.updateCfg(obj.cfg,varargin{:});
            obj.t0stamp = obj.isoTimestamp();
            obj.setGlobalStyles();
        end

        function hFig = newFigure(obj,name)
            obj.idx = obj.idx + 1;
            figName = name;

            hFig = figure( ...
                'Name', figName, ...
                'NumberTitle', 'off', ...
                'Visible', obj.bool2onoff(obj.cfg.visible), ...
                'Units', 'normalized', ...
                'OuterPosition',[0 0 1 1]);

            plotbrowser('on');
            ax = axes('Parent',hFig); hold(ax,'on'); grid(ax,'on');
            title(ax,name);

            obj.handles.(figName) = hFig;
        end

        function exportAll(obj)
            set(findall(0,'Type','figure'),'Visible','on');
            if (~ obj.cfg.export_figs)
                return
            end

            if ~isfolder(obj.cfg.exportDir), mkdir(obj.cfg.exportDir); end
%             obj.writeMetadata();

            flds = fieldnames(obj.handles);

            pause(1)

            for k = 1:numel(flds)
                h = obj.handles.(flds{k});
                plotedit(h,'off')
                
%                 tag = sprintf('%s_%s',flds{k},obj.t0stamp);
                tag = flds{k};

                for fmt = obj.cfg.formats
                    fn = fullfile(obj.cfg.exportDir,tag+"."+fmt);
                    opts = {'Resolution',obj.cfg.dpi};
                    if obj.cfg.transparentBG, opts = [opts,{'BackgroundColor','none'}]; end
                    exportgraphics(h,fn,opts{:});
                end
                if obj.cfg.saveFig
                    savefig(h,fullfile(obj.cfg.exportDir,tag+".fig"));
                end
            end
        end

        %% Convenience ----------------------------------------------------
        function setCurrent(obj,fig)
            if ischar(fig) || isstring(fig)
                figure(obj.handles.(char(fig)));
            elseif ishghandle(fig,'figure')
                figure(fig);
            else
                error('Input must be a figure handle or stored name.');
            end
        end

        function names = listFigures(obj)
            names = fieldnames(obj.handles);
        end

        function closeFigure(obj,fig)
            if ischar(fig) || isstring(fig), fig = char(fig); end
            if ischar(fig)
                if isfield(obj.handles,fig)
                    close(obj.handles.(fig));
                    obj.handles = rmfield(obj.handles,fig);
                end
            elseif ishghandle(fig,'figure')
                close(fig);
                flds = fieldnames(obj.handles);
                for k = 1:numel(flds)
                    if obj.handles.(flds{k}) == fig
                        obj.handles = rmfield(obj.handles,flds{k}); break;
                    end
                end
            end
        end

        function closeAll(obj)
            flds = fieldnames(obj.handles);
            for k = 1:numel(flds), close(obj.handles.(flds{k})); end
            obj.handles = struct(); obj.idx = 0;
        end
    end

    %% PRIVATE ============================================================
    methods (Access=private)
        function cfg = defaultCfg(~)
            cfg = struct( ...
                'exportDir', fullfile(pwd,'exports'), ...
                'visible', true, ...
                'export_figs', false, ...
                'formats', "png", ...
                'saveFig', false, ...
                'dpi', 300, ...
                'transparentBG', false, ...
                'colourOrder', lines(7), ...
                'lineStyleOrder','-' ...
                );
        end

        function cfg = updateCfg(~,cfg,varargin)
            if mod(numel(varargin),2)~=0, error('Name-Value pairs required'); end
            for v = 1:2:numel(varargin)
                key = varargin{v}; val = varargin{v+1};
                if isfield(cfg,key), cfg.(key) = val; end
            end
        end

        function setGlobalStyles(obj)
            set(groot,'defaultAxesColorOrder',obj.cfg.colourOrder);
            set(groot,'defaultAxesLineStyleOrder',strjoin({'-','--','.-',':'},'|'));
        end

        function writeMetadata(obj)
            fn = fullfile(obj.cfg.exportDir,"config_"+obj.t0stamp+".json");
            fid = fopen(fn,'w');
            if fid~=-1
                fprintf(fid,'%s',jsonencode(obj.cfg,'PrettyPrint',true));
                fclose(fid);
            end
        end
    end

    %% STATIC PRIVATE =====================================================
    methods (Static, Access=private)
        function s = bool2onoff(b), s = ["off","on"]; s = s(b+1); end
        function ts = isoTimestamp, ts = datestr(now,'yyyymmdd_HHMMSS'); end
    end
end