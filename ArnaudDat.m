

classdef ArnaudDat < handle
    
    properties
        % Raw data & related metadata
        dat             % Original raw data; last 2 columns are the trigger channels
        t
        dt
        do_part
        do_MUA
        
        
        % Data blocked into trials & related metadata
        dat_tr          % Data trial blocks
        t_tr            % Time trial blocks
        is_AP           % Boolean, true/false for AP pulse
        Ntrials         % Number of trials composing this data.
        APname = ''     % Text information about whether AP data or not
        
        % Figure metadata
        visible = 'on';        % 'on' or 'off' for setting new figures to be visible
        fname = [];
        gcf = {}                   % Array of figure handles
        figtype = {}               % Array of figure types
        myylabel;
        
        % Other metadata
        pp1
        pp2
        trig1
        trig2
        trigchan1
        trigchan2
        ind
        Nchan
        
        
    end
    
    properties (Access = private)
        
    end
    
    methods
        function s = ArnaudDat(d,fname,do_MUA,do_part,visible)
            % Constructor!
            
            % Set up defaults
            if nargin < 5
                visible = [];
            end
            if nargin < 4
                do_part = [];
            end
            if nargin < 3
                do_MUA = [];
            end
            if nargin < 2
                fname = [];
            end
            if isempty(visible)
                visible = 'on';
            end
            if isempty(do_part)
                do_part = true;
                do_part = false;
            end
            if isempty(do_MUA)
                do_MUA = false;
            end
            
            
            if do_part
                ind = 1:10000;               % Do part of dataset (~4 cycles)
                % ind = 1:30000;             % Do part of dataset (~12 cycles of pulse train)
            else
                if do_MUA; ind = 1:length(d.cntm); % Do full lenght of dataset              
                else; ind = 1:length(d.cnte); % Do full lenght of dataset             
                end
            end
            
            % Pull in data
            vars_pull(d);
            
            % Time variables
            dt = mode(diff(time));
            t = time(ind);
            
            
            % Extract pulse train times
            pp1 = anatrig{1};
            pp2 = anatrig{2};
            pp1 = pp1(pp1 < max(ind));
            pp2 = pp2(pp2 < max(ind));
            
            % Convert to triggers marking the start of each ppstim train
            thresh = 1000;
            trig1 = pp2trig (pp1,thresh);
            trig2 = pp2trig (pp2,thresh);
            
            % Format the data we're looking at
            if do_MUA
                dat0 = cntm(:,ind)';
            else
                dat0 = cnte(:,ind)';
            end
            trigchan1 = dat0(:,trigchan(1));
            trigchan2 = dat0(:,trigchan(2));
            dat = dat0(:,chansets(1):chansets(2));
            
            % dat = diff(dat,2,2);       % 2nd difference to approx laplacean
            dat = cat(2,dat,trigchan1,trigchan2);   % Add back in pulse channels
            Nchan = size(dat,2);
            
            % Pack in all data
            s.dat = dat;
            s.t = t;
            s.dt = dt;
            s.pp1 = pp1;
            s.pp2 = pp2;
            s.trig1 = trig1;
            s.trig2 = trig2;
            s.trigchan1 = trigchan1;
            s.trigchan2 = trigchan2;
            s.ind = ind;
            s.Nchan = Nchan;
            [~,fname,~] = fileparts(fname);
            s.visible = visible;
            
            if do_MUA           % Append MUA information
                fname = [fname '_MUA'];
            else
                fname = [fname '_CSD'];
            end
            if do_part          % Append part information
                fname = [fname '_part'];
            else
                fname = [fname '_full'];
            end
            s.fname = fname;
            s.do_part = do_part;
            if do_MUA
                s.myylabel = 'MUA firing rate';
            else
                s.myylabel = 'CSD amp';
            end
        end
        
        function s = plot_raw(s,do_shift)
            
            if nargin < 2
                do_shift = false;
            end
            
            if ~s.do_part
                error('Data set too large to plot. plot_raw only runs with do_part = false.');
            end
            
            d = s.dat;

            [~,~,method_name] = fileparts(calledby(0)); method_name = method_name(2:end);
            
            if do_shift
                s.gcf{end+1} = figw('visible',s.visible); 
                s.figtype{end+1} =  [method_name '_shift'];
                plot_with_shift(s.t,d,s);
                
            else
                show_trigger_traces = false;
                % Plot data
                s.gcf{end+1} = figw('visible',s.visible); 
                %plot_default(s.t,d(:,[1:3,end-1,end]),s);
                
                plot_default(s.t,d(:,1:end),s,show_trigger_traces);
                s.figtype{end+1} = [method_name '_noshift'];          % Cheesy way of getting current method. Maybe there's a specific command for this...
                
                
                if show_trigger_traces
                    % Add tick markers
                    hold on;
                    % plot(t(pp1),dat(pp1,:),'k.');
                    % plot(t(pp2),dat(pp2,:),'r.');
                    plot(s.t(s.trig1),d(s.trig1,s.Nchan-1:s.Nchan),'k.');
                    plot(s.t(s.trig2),d(s.trig2,s.Nchan-1:s.Nchan),'r.');
                end
            end
        end
        
        function s = overlay_raw(s,channels_to_plot,upscale_val,varargin)
            
            
            if nargin < 3
                upscale_val = [];
            end
            if nargin < 2
                channels_to_plot = [];
            end

            
            plot_overlay(s.t,s.dat,s,channels_to_plot,upscale_val,varargin{:});
            
        end
        
        function s = imagesc_raw(s)
            
            time = s.t;
            d = s.dat;
            
            if ~s.do_part
                error('Data set too large to plot. plot_raw only runs with do_part = false.');
            end
            
            
            % Find y limits.
            yl = [min(min(d(:,1:s.Nchan-2))),max(max(d(:,1:s.Nchan-2)))];       % Subtract 2 for the two spiking channels
            
            
            % Plot Imagesc
            s.gcf{end+1} = figw('visible',s.visible);  imagesc([time(1), time(end)],[],d'); axis xy; caxis([yl(1) yl(2)]); %colorbar
            
            % Add tic labels
            addTicks(time,d(:,s.Nchan-1),'r');
            addTicks(time,d(:,s.Nchan),'b');
    
            [~,~,method_name] = fileparts(calledby(0));  method_name = method_name(2:end);
            s.figtype{end+1} = method_name;
            
            % Add title
            xlabel('time (ms)');
            ylabel('channel #');
            title (format_title([s.fname ' # ch=' num2str(s.Nchan - 2) ' # tr=' num2str(s.Ntrials)]));
            
            colormap jet
        end
        
        function s = calc_trialblocks(s,do_AP)
            % Reshape the data into trials
            
            if nargin < 2
                do_AP = [];
            end
            
            if isempty(do_AP); do_AP = true; end
            
            if do_AP
                s.is_AP = true;
                mytrig = s.trig1;
                s.APname = ['_AP'];          % Aperiodic pulse
            else
                s.is_AP = false;
                mytrig = s.trig2;
                s.APname = ['_PP'];          % Periodic pulse
            end
            
            d = s.dat;
            
            tr_start = round(-100/s.dt);         % trial starts at -100ms
            tr_end = round(600/s.dt - 1);        % Ends at 600 ms
            
            trind = repmat([tr_start:tr_end]',1,length(mytrig));         % Block of indices for selecting trials
            trind = trind + repmat(mytrig(:)',size(trind,1),1);
            
            % Remove any trind columns that exceeed the size of dat2
            trind = trind(:, max(trind) <= size(d,1));
            
            % dch1 = dat(:,1);
            % t1 = dch1(trind);
            
            dat2 = permute(d,[1,3,2]);
            d3 = dat2(trind,:,:);
            sz = size(trind);
            d3 = reshape(d3,[sz(1), sz(2), s.Nchan]);
            
            
            s.dat_tr = d3;
            s.t_tr = [tr_start:tr_end] * s.dt;
            
            s.Ntrials = size(d3,2);
            
            
        end
        
        function s = plot_travg(s,do_shift)

            
            if nargin < 2
                do_shift = false;
            end
            dtr = s.dat_tr;
            dtr = squeeze(mean(dtr,2));
            
            s.gcf{end+1} = figw('visible',s.visible); 
            [~,~,method_name] = fileparts(calledby(0)); method_name = method_name(2:end);
            if do_shift
                if s.Ntrials > 0; plot_with_shift(s.t_tr,dtr,s); end
                s.figtype{end+1} = [method_name '_shift'];
            else 
                if s.Ntrials > 0; plot_default(s.t_tr,dtr,s); end
                s.figtype{end+1} = [method_name '_noshift'];
            end
        end
        
        function s = overlay_travg(s,channels_to_plot,upscale_val,varargin)
            
            if nargin < 3
                upscale_val = [];
            end
            if nargin < 2
                channels_to_plot = [];
            end

            dtr = s.dat_tr;
            dtr = squeeze(mean(dtr,2));
            

            [~,~,method_name] = fileparts(calledby(0)); method_name = method_name(2:end);
            
            plot_overlay(s.t_tr,dtr,s,channels_to_plot,upscale_val,varargin{:});
            s.figtype{end+1} = [method_name '_overlay'];
        end
        
        function s = imagesc_travg(s)
            
            
            
            dtr = s.dat_tr;
            dtr = squeeze(mean(dtr,2));
            ttr = s.t_tr;
            
            if s.Ntrials > 0
                yl = [min(min(dtr(:,1:s.Nchan-2))),max(max(dtr(:,1:s.Nchan-2)))];       % Subtract 2 for the two spiking channels
                s.gcf{end+1} = figw('visible',s.visible);  imagesc([ttr(1) ttr(end)],[],dtr'); caxis([yl(1) yl(2)]); %colorbar
            else
                s.gcf{end+1} = figw('visible',s.visible);
            end
            [~,~,method_name] = fileparts(calledby(0)); method_name = method_name(2:end);
            s.figtype{end+1} = method_name;
            axis xy
            xlabel('time (ms)');
            ylabel('channel #');
            title (format_title([s.fname ' # ch=' num2str(s.Nchan - 2) ' # tr=' num2str(s.Ntrials)]));
            
            colormap jet
            
        end
        
        function save_openFigs(s,mypath)
            if nargin < 2
                mypath = [];
            end
            if isempty(mypath)
                mypath = '.';
            end
            for i = 1:length(s.gcf)
                curr_handle = s.gcf{i};
                if ishandle(curr_handle)
                    filename = [fullfile(mypath,[s.fname '_' s.APname '_' s.figtype{i}]) '.png'];
                    if exist([filename],'file'); error('File %s already exists! Exiting...',filename); end
                    set(curr_handle,'PaperPositionMode','auto');
                    fprintf('Saving %s ... ',filename);
                    tic; print(curr_handle,'-dpng','-r75','-opengl',filename);toc
                    
                end
            end
        end
        
        function clearFigs(s)
            close all;
            s.gcf = {};
            s.figtype = {};
        end
        
    end
    
end


function addTicks (t,d,colour)

    threshold = 2;
    lw = 20;
    
    ind = d > threshold;
    d(ind) = 1;
    d(~ind) = nan;
    
    yl = ylim;
    
    hold on; plot(t,d*yl(2)-0.01,[colour '--.'],'LineWidth',lw);
    

end

function plot_default(t,d,s,show_trigger_traces)

    if nargin < 4
        show_trigger_traces = [];
    end
    
    if isempty(show_trigger_traces)
        show_trigger_traces = true;
    end
    
    if show_trigger_traces
        plot(t,d);
    else
        plot(t,d(:,1:end-2));
    end
    
    xlabel('time (ms)');
    %set(gca,'YTick',[]);
    
    ylabel(s.myylabel);
    xlabel('time (ms)');
    title (format_title([s.fname ' # ch=' num2str(s.Nchan - 2) ' # tr=' num2str(s.Ntrials)]));
    
    addTicks(t,d(:,s.Nchan-1),'r');
    addTicks(t,d(:,s.Nchan),'b');

end

function plot_overlay(t,d,s, channels_to_plot,upscale_val,varargin)

    if nargin < 5
        upscale_val = [];           % Amplitude scaling of MUA traces.
                                    % 1 = MUA trace 1's max equals MUA
                                    % trace'e 2's min. Larger numbers
                                    % provide greater overlap
    end
    
    if nargin < 4
        channels_to_plot = [];
    end
    
    if isempty(upscale_val)
        upscale_val = 1;
    end
    if isempty(channels_to_plot)
        channels_to_plot = 1:s.Nchan;
    end

    % Adjust amplitude and vertical spacing of data
    d = ampscore(d);                                                % Adjust data so that the minimum is at zero and the maximum is at 1
    d = ampscale(d,min(diff(channels_to_plot))*upscale_val);        % Scale amplitude to be equal to the smallest spacing between chosen channels
    
    % Shift each channel by the channel number (e.g. so channel 1 will
    % start at 1, channel 2 at 2, etc.)
    shifts = 1:s.Nchan;
    d = d + repmat(shifts,size(d,1),1);
    
    hold on;
    plot(t,d(:,channels_to_plot),varargin{:});
    
    %xlabel('time (ms)');
    %set(gca,'YTick',[]);
    
    %ylabel(s.myylabel);
    %xlabel('time (ms)');
    %title (format_title([s.fname ' # ch=' num2str(s.Nchan - 2) ' # tr=' num2str(s.Ntrials)]));

end
function plot_with_shift(t,d,s)
    Nchan = s.Nchan;
    plott_matrix3D(t,d(:,1:Nchan-2),'active_dim',3,'do_shift',mean(std(d,[],1))*2);
    hold on;
    addTicks(t,d(:,Nchan-1),'r');
    addTicks(t,d(:,Nchan),'b');
    
    xlabel('time (ms)');
    set(gca,'YTick',[]);
    
    xlabel('time (ms)');
    title (format_title([s.fname ' # ch=' num2str(s.Nchan - 2) ' # tr=' num2str(s.Ntrials)]));
    
end

function varargout= figw(varargin)
    varargout = cell(1,nargout);
    [varargout{1:nargout}] = figure('Color','w','Position',[41 391 1400 407],varargin{:});

end

function str = format_title(str)
    str = strrep(str,'_',' ' );
end


function X = ampscore(X)
    % Adjust data so that the minimum is at zero and the maximum is at 1
    sz = size(X);
    for i = 1:sz(2)
        Xc = X(:,i);
        Xc = Xc ./ (max(Xc) - min(Xc));     % Normalize to 1
        Xc = Xc - min(Xc);                  % Starts at zero.
        X(:,i) = Xc;
    end
end


function X = ampscale(X,val)
    % Adjust data so that the minimum is at zero and the maximum is at 1
    sz = size(X);
    for i = 1:sz(2)
        Xc = X(:,i);
        Xc = Xc * val;      % Scale amplitude of this column to value
        X(:,i) = Xc;
    end
end