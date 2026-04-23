%% ENF Signature Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:
% This script implements the signal processing pipeline for the application
% "Geolocation for Semiconductor Chips using Electric Network Frequency Signatures".
%
% The primary goal is to validate the feasibility of extracting Electric
% Network Frequency (ENF) signatures from DC-powered hardware by comparing
% a sensed trace against a ground-truth reference.
%
% Analysis Steps:
%   1.  INPUTS: The script takes two time-synchronized inputs (both must be .wav):
%       - A 'ground-truth reference trace' captured from the AC mains. (.wav file format)
%       - A 'sensed trace' captured from the experimental FPGA board (.wav file format).
%         This can be either an ambient EM trace or a board power trace.
%   
%   2.  SPECTROGRAM GENERATION: The Short-Time Fourier Transform (STFT) is
%       applied to both traces to generate high-resolution spectrograms,
%       visualizing the harmonic content over time.
%   
%   3.  ENF ESTIMATION: A weighted average PMF is used to extract the instantaneous ENF
%       signature from each spectrogram.
%   
%   4.  OUTPUTS & CORRELATION: Finally, the script compares the ENF signature
%       from the sensed trace against the signature from the ground-truth
%       reference. It calculates the Pearson correlation coefficient and
%       generates the final temporally aligned plots to visually and quantitatively
%       assess the match.
%
% Dependencies:
%   proc_enf_analysis - Pre-compiled ENF extraction and correlation function.
%                       Must be on the MATLAB path before running this script.
%
% Usage:
%   Run from the exp_scripts/ directory. The artifact data root is resolved
%   automatically relative to this script file's location.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Script Configuration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ENF Analysis Configuration Flags:
% The parameters in this section control what is analyzed, how the STFT spectrogram
% is computed, how the ENF frequency is estimated, and how results are displayed.
% Update the four variables in the 'Input trace file information' block for each
% new run; all other parameters can typically remain at their default values.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input trace file information
% This script is designed to analyze ONE pair of traces at a time. Manually
% update the variables below for EACH run. If both traces are in the same
% leaf folder, set file_1_path and file_2_path to the same value.
%
% 1. `file_1_path`: Relative path ending in "/" for the AC mains reference.
% 2. `file_2_path`: Relative path ending in "/" for the sensed trace.
% 3. `file_1_name`: AC mains reference file name without extension.
% 4. `file_2_name`: Sensed trace file name without extension.
%
% --- EXAMPLES ---
%
% % Example 1: Analyzing a 'TREND' experiment trace (WK01/WEND/SUN/EVEN/T01)
% file_1_path = "../exp_inputs/TREND/WK01/WEND/SUN/EVEN/T01/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 2: Analyzing an FPGA ambient EM trace (FPGA/EM_TRACES/CW305/DC_1)
% file_1_path = "../exp_inputs/FPGA/EM_TRACES/CW305/DC_1/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";
%
% % Example 3: Analyzing an FPGA board power trace (FPGA/POW_TRACES/CW305/DC_1)
% file_1_path = "../exp_inputs/FPGA/POW_TRACES/CW305/DC_1/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_pow_trace_dc_a";  % or "fpga_pow_trace_dc_b"
%
% % Example 4: Analyzing a 'MULTI' (US_60) experiment trace (MULTI/US_60/AUG/WED/T02)
% % File names differ from the standard naming convention in this folder.
% file_1_path = "../exp_inputs/MULTI/US_60/AUG/WED/T02/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac_egrid_citya_lab";
% file_2_name = "fpga_em_trace_dc_egrid_citya_lab";
%
% % Example 5: Analyzing a 'MULTI' (DE_50) experiment trace (MULTI/DE_50/AUG/WED/T01)
% file_1_path = "../exp_inputs/MULTI/DE_50/AUG/WED/T01/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac_citya_lab";
% file_2_name = "fpga_em_trace_dc_citya_lab";
%
% % Example 6: Analyzing an 'SRV_L' experiment trace (SRV_L/SPSU/WED/T05)
% file_1_path = "../exp_inputs/SRV_L/SPSU/WED/T05/";
% file_2_path = file_1_path;
% file_1_name = "mains_pow_trace_ac";
% file_2_name = "fpga_em_trace_dc";

% Set the variables below to point to the trace pair you want to analyze.
% file_1_path = "../exp_inputs/FPGA/EM_TRACES/CW305/DC_1/";
% file_2_path = file_1_path;
file_1_path = "../exp_inputs/MULTI/US_60/OCT/WED/T01/";
file_2_path = file_1_path;
file_1_name = "mains_pow_trace_ac_egrid_citya_lab";
file_2_name = "fpga_em_trace_dc_egrid_citya_lab";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spectrogram Computation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fundamental power grid frequency settings
nominal_freq_arr      = [50 60];

% --- Multi-Location Experiment Note (DE_50 / US_60) ---
% When analyzing traces from the 'MULTI' directory, you MUST manually
% update the `nominal_freq_*` variables below to match the fundamental
% frequency of the grid from where the trace was recorded.
% - For US_60 (United States) traces: set nominal_freq = nominal_freq_arr(2); (60 Hz)
% - For DE_50 (Germany) traces:     set nominal_freq = nominal_freq_arr(1); (50 Hz)

%For reference trace
nominal_freq_1        = nominal_freq_arr(2);
harmonics_arr_1       = (1:7)*nominal_freq_1;

%For sensed trace
nominal_freq_2        = nominal_freq_arr(2);
harmonics_arr_2       = (1:7)*nominal_freq_2;

% STFT compute param settings
frame_size_arr      = (1:12)*1000;
frame_size          = frame_size_arr(8);                %8000ms window
nfft_arr            = 2.^(10:20);
nfft                = nfft_arr(6);                      %2^15 = 32768 pts
overlap_size_arr    = 0:0.1:0.9;
overlap_size        = overlap_size_arr(1)*frame_size;   %non-overlapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency Estimation Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Frequency to be estimated
trace_1_est_freq = harmonics_arr_1(1);                    %1st harmonic= 60Hz
trace_2_est_freq = harmonics_arr_2(1);                    %1st harmonic= 60Hz

%Frequency estimation method options:
%1: weighted average (pmf power = 3)
%2: spectrum combining (quad interp)
trace_1_freq_est_method = 1;   %For reference trace
trace_1_freq_est_spec_comb_harmonics = [60 120 180 240 300 360 420];

trace_2_freq_est_method = 1;   %For sensed trace
trace_2_freq_est_spec_comb_harmonics = [60 120 180 240 300 360 420];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MATLAB Plotting Parameters:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The colormap for the figures is set to jet.
set(0,'DefaultFigureColormap', jet)

%Titles for MATLAB plots
trace_1_plot_title = "Reference AC Mains Power Trace";
trace_2_plot_title = "Sensed FPGA Trace";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% File Path Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build the full .wav file paths for the reference and sensed traces by
% concatenating each configured leaf-folder path, file base name, and the
% '.wav' extension. Both inputs must already be in .wav format; no format
% conversion is performed here.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
full_file_1_path_wav = file_1_path + file_1_name + ".wav";
full_file_2_path_wav = file_2_path + file_2_name + ".wav";
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Main Analysis: ENF Extraction and Correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Runs the full ENF analysis pipeline on the configured reference and sensed
% trace pair. All signal processing is handled by proc_enf_analysis, which:
%   - Computes the STFT spectrogram for both traces
%   - Extracts the instantaneous ENF time series from each spectrogram
%   - Temporally aligns the two ENF traces and computes Pearson correlation
%   - Generates all intermediate spectrograms and the final correlation plot
% The returned value is the best (temporally aligned) Pearson correlation (r).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\nINFO: ENF Analysis on complete reference and sensed traces ...\n');
matched_corr_val_full = proc_enf_analysis(full_file_1_path_wav, full_file_2_path_wav, ...
                                          nfft, frame_size, overlap_size, ...
                                          harmonics_arr_1, nominal_freq_1, harmonics_arr_2, nominal_freq_2,...
                                          trace_1_freq_est_method, trace_1_est_freq, trace_1_freq_est_spec_comb_harmonics, ...
                                          trace_2_freq_est_method, trace_2_est_freq, trace_2_freq_est_spec_comb_harmonics, ...
                                          trace_1_plot_title, trace_2_plot_title, true);

fprintf('\nINFO: Pearsons Correlation Coefficient for the temporally aligned reference and sensed traces (complete trace): %.8f \n', matched_corr_val_full);
fprintf('\nINFO: ENF Analysis complete.\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
