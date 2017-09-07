


%% Description
% This is the script for plotting basic CSD and MUA data collected by
% Arnaud. Has some more advanced options than just myscript. Uses class.
% Generally, most of my work is done using this script now. 

%% Setup path info

startup

mypath = fullfile('..','Shared_Volume_1');


%% Load a single
% Load some data
% fname = 'WI03021@o.mat';
% fname = 'WI03022@o.mat';
fname = 'WI10030@o.mat';
d = load(fullfile(mypath,fname));

clearvars -except d mypath

vars_pull(d)


%% Plot a single CSD file

clear ad
do_MUA = false;
do_part = false;
do_AP = true;

ad = ArnaudDat(d,'WI03021@o.mat',do_MUA,do_part);
ad.calc_trialblocks(do_AP);


% ad.plot_raw(true);
% ad.imagesc_raw;


% ad.plot_travg(true);
ad.imagesc_travg;


%% Overlay MUA 
do_MUA = true;

admua = ArnaudDat(d,'WI03021@o.mat',do_MUA,do_part);
admua.calc_trialblocks(do_AP);

average_unplotted_traces = true;

% Plot trial averages
ad.imagesc_travg; admua.overlay_travg(1:3:21,2,average_unplotted_traces,'k');

% Plot raw
% ad.imagesc_raw; admua.overlay_raw(1:3:21,1,average_unplotted_traces,'k');
% ad.imagesc_raw; ad.overlay_raw(1:2:21,2,average_unplotted_traces,'k');


%% Figs1: Load and save all files! (No overlays; loopy version)

n = dir(fullfile(mypath,'*.mat'));
for i = 1:length(n)
    
    currfile = n(i).name;
    fprintf('Loading file %s...', currfile);
    d = load(fullfile(mypath,currfile));
    fprintf('Done.\n');
    for do_MUA = [false, true]
        for do_part = [false, true]
            clear ad
            ad = ArnaudDat(d,currfile,do_MUA,do_part);
            ad.visible = 'off';
            
            if do_part
                do_shift_temp = true; ad.plot_raw(do_shift_temp);
                do_shift_temp = false; ad.plot_raw(do_shift_temp);
                %ad.imagesc_raw;
                
                ad.save_openFigs;
                ad.clearFigs;
                
            else
                for do_AP = [false, true]

                    ad.calc_trialblocks(do_AP);

                    %ad.plot_travg(do_shift);
                    ad.imagesc_travg;
                    
                    ad.save_openFigs;
                    ad.clearFigs;

                end
            end
        end
    end
    
end



%% Figs2: Load and save all files! (With MUA overlaid on CSD)

n = dir(fullfile(mypath,'*.mat'));

do_shift = true;
average_unplotted_traces = true;

for i = 1:length(n)
    
    currfile = n(i).name;
    fprintf('Loading file %s...', currfile);
    d = load(fullfile(mypath,currfile));
    fprintf('Done.\n');
    
    % First, load part of data
%     clear ad admua
%     do_part = true;
%     ad = ArnaudDat(d,currfile,false,do_part,'off');
%     admua = ArnaudDat(d,currfile,true,do_part,'off');
%     
%     ad.plot_raw(false);
%     ad.imagesc_raw; admua.overlay_raw(1:3:21,1,average_unplotted_traces,'k');
%     ad.imagesc_raw; ad.overlay_raw(1:2:21,2,average_unplotted_traces,'k'); ad.figtype{end} = [ad.figtype{end} '2'];  % Rename figure so don't get duplicates
%     ad.save_openFigs;
%     ad.clearFigs;
    
    % Next, plot all trial averaged data
    clear ad admua
    do_part = false;
    for do_AP = [false, true]
        ad = ArnaudDat(d,currfile,false,do_part,'off');
        admua = ArnaudDat(d,currfile,true,do_part,'off');
        
        ad.calc_trialblocks(do_AP);
        admua.calc_trialblocks(do_AP);

        ad.imagesc_travg; admua.overlay_travg(1:3:21,2,average_unplotted_traces,'k');
        
        ad.save_openFigs;
        ad.clearFigs;
        
    end
end


