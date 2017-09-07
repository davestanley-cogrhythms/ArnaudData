
%% Description
% This is the script for plotting basic CSD and MUA data collected by
% Arnaud.

%% 
startup
mypath = fullfile('..','Shared_Volume_1');

% Load some data
d = load(fullfile(mypath,'WI03021@o.mat'));
% d = load(fullfile(mypath,'WI03022@o.mat'));
clearvars -except d

vars_pull(d)

plot_raw_on = true;

%% PSD / Spectrogram plots

figure; plott_spect(t/1000,dat(:,10),'axis_lims',[0 50],'logplot',1,'Nwind',2000)
figure; plott_all(t/1000,dat(:,10),'axis_lims',[0 50],'logplot',1,'psd_on',1,'fsubplot',@subplotrows,'Nwind',2000)

%% do some basic plots
% close all
ind = 1:90000;               % Do part of dataset
% ind = 1:length(cntm);               % Do full lenght of dataset
dt = mode(diff(time));

% Format the data we're looking at
dat0 = cnte(:,ind)';
trigchan1 = dat0(:,trigchan(1));
trigchan2 = dat0(:,trigchan(2));

dat = dat0(:,chansets(1):chansets(2));

% dat = diff(dat,2,2);       % 2nd difference to approx laplacean
dat = cat(2,dat,trigchan1,trigchan2);   % Add back in pulse channels
Nchan = size(dat,2);


t = time(ind);
yl = [min(min(dat(:,1:Nchan-2))),max(max(dat(:,1:Nchan-2)))];       % Subtract 2 for the two spiking channels

% Extract pulse train times
pp1 = anatrig{1};
pp2 = anatrig{2};
pp1 = pp1(pp1 < max(ind));
pp2 = pp2(pp2 < max(ind));

% Convert to triggers marking the start of each ppstim train
thresh = 1000;
trig1 = pp2trig (pp1,thresh);
trig2 = pp2trig (pp2,thresh);

% Plot data
if plot_raw_on
    figure; plot(t,dat(:,[1:3,end-1,end]));
end

% Add on PP time markers info

if plot_raw_on
    hold on;
    % plot(t(pp1),dat(pp1,:),'k.');
    % plot(t(pp2),dat(pp2,:),'r.');
    plot(t(trig1),dat(trig1,:),'k.');
    plot(t(trig2),dat(trig2,:),'r.');
end

%% Plot Imagesc
% figure; imagesc([time(1), time(end)],[],dat'); axis xy; caxis([yl(1) yl(2)]); colorbar


%% Reshape the data into trials

tr_start = round(-100/dt);         % 600 ms trial lengths
tr_end = round(600/dt);         % 600 ms trial lengths

trind = repmat([tr_start:tr_end]',1,length(trig1));         % Block of indices for selecting trials
trind = trind + repmat(trig1(:)',size(trind,1),1);

% Remove any trind columns that exceeed the size of dat2
trind = trind(:, max(trind) <= size(dat,1));

% dch1 = dat(:,1);
% t1 = dch1(trind);

dat2 = permute(dat,[1,3,2]);
dat3 = dat2(trind,:,:);
sz = size(trind);
dat3 = reshape(dat3,[sz(1), sz(2), Nchan]);


dat3 = squeeze(mean(dat3,2));

yl = [min(min(dat3(:,1:Nchan-2))),max(max(dat3(:,1:Nchan-2)))];       % Subtract 2 for the two spiking channels
figure; imagesc(dat3'); caxis([yl(1) yl(2)]); colorbar
axis xy
title (['# trials = ' num2str(size(dat3,2))]);

