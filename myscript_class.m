


%% Description
% This is the script for plotting basic CSD and MUA data collected by
% Arnaud. Has some more advanced options than just myscript. Uses class.
% Generally, most of my work is done using this script now. 

%% Setup path info

startup

mypath = fullfile('..','Shared_Volume_1');


%% Load a single
% Load some data
d = load(fullfile(mypath,'WI03021@o.mat'));
% d = load(fullfile(mypath,'WI03022@o.mat'));
clearvars -except d

vars_pull(d)


%% Plot a single file

clear ad
do_MUA = false;
do_part = true;
do_AP = false;

ad = ArnaudDat(d,'WI03021@o.mat',do_MUA,do_part);
ad.calc_trialblocks(do_AP);


ad.plot_raw(true);
ad.imagesc_raw;


ad.plot_travg(true);
ad.imagesc_travg;


%% Load and save all files!

n = dir(fullfile(mypath,'*.mat'));

do_shift = true; 

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

