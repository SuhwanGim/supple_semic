%% DAY1
clc;
clear;
close all;
%% SETUP: PARAMETER
%ip = '115.145.189.133';
ip = '192.168.0.3'; port = 20121; addpath(genpath(pwd));

%% EXPERIMENT %%
%% 1. Calibration
% 'SEMP000~'
SID = 입력하세요 %Participants
calibration(SID, ip, port,'joystick'); %'test'
%% SETUP: Load calibration data & SkinSite sequences
reg = load_cali_results(SID); disp(reg.skinSite_rs);
%% Main STUDY % thermode_test(('biopac1','fmri') ... ... )
%% 1. RUN1
thermode_test(SID, 1, ip, port, reg,'biopac1','joystick','test');
%% 2. RUN2
thermode_test(SID, 2, ip, port, reg,'biopac1','joystick','test');
%% 3. RUN3
thermode_test(SID, 3, ip, port, reg,'biopac1','joystick','test');
%% 4. RUN4
thermode_test(SID, 4, ip, port, reg,'biopac1','joystick','test');




%% SEND DATA
SEMIC_Senddata(SID,'mri')