


clear

addpath(genpath('G:\My Drive\GitHub\fooof_mat\'))
% pyenv('ExecutionMode', 'OutOfProcess');

% pyenv('Version', 'C:\Users\dtf8829\AppData\Local\anaconda3\python.exe')

data = load('H:\SUMDAT\biomarkCon_7_7.mat').data; 
frex = logspace(log10(2),log10(80),100);

psd = squeeze(mean(data.power(32,:,:),3));

% FOOOF settings
settings = struct();  % Use defaults
f_range = [1, 30];

% Run FOOOF
fooof_results = fooof(frex, psd, f_range, settings);



% % Define settings for FOOOF
settings.peak_width_limits = [1, 6];
settings.max_n_peaks = 6;
settings.min_peak_height = 0.1;
settings.peak_threshold = 2;
% 
% % Example frequency and PSD data
% frex = linspace(1, 40, 100);  % Frequency range from 1 to 40 Hz
% psd = rand(1, 100);           % Random example PSD data
% f_range = [3, 30];            % Frequency range to fit FOOOF model

% Import the fooof module
py_fooof = py.importlib.import_module('fooof');
np = py.importlib.import_module('numpy');  % Import NumPy

% Convert MATLAB arrays to NumPy arrays
frex_np = np.array(frex);      % Convert frequencies to NumPy array
psd_np = np.array(psd);        % Convert PSD data to NumPy array
f_range_np = np.array(f_range); % Convert frequency range to NumPy array

% Create an instance of the FOOOF model with settings
fm = py_fooof.FOOOF();  % Using default settings for simplicity

% Fit the model to the data
fooof_results = fm.fit(frex_np, psd_np, f_range_np);



sumdat = dir('H:\SUMDAT\');
sumdat(1:2,:) = []; 


frex = logspace(log10(2),log10(80),100);
for ii = 1:length(sumdat)
   
    ii
    data = load([sumdat(ii).folder '/' sumdat(ii).name]).data; 

    psd = squeeze(mean(data.power(32,:,:),3));

    T = table(psd', frex');

    % Optionally, you can name the columns for clarity
    T.Properties.VariableNames = {'power', 'frex'};
    T.slope = ones(100,1)*data.slopeValsLog(32,1);
    T.offset = ones(100,1)*data.slopeValsLog(32,2);
    
    % Save the table as a CSV file
    writetable(T, [data.key '.csv']);
   
end

















