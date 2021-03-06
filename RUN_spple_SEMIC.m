%% DAY1
clc;
clear;
close all;
%% SETUP: PARAMETER
%ip = '115.145.189.133';
ip = '192.168.0.3'; port = 20121; addpath(genpath(pwd));
type = rnd_semic_type;


%% EXPERIMENT %%
%% 1. Calibration
% 'SEMP000~'
SID = 'SEMP002' %Participants % pSEM001~
calibration(SID, ip, port,'joystick'); %'test'
%% SETUP: Load calibration data & SkinSite sequences
reg = load_cali_results(SID); disp(reg.skinSite_rs);
%% Main STUDY % thermode_test(('biopac1','fmri') ... ... )
%% 1. RUN1r
thermode_test(type, SID, 1, ip, port, reg,'joystick');
%% 2. RUN2
thermode_test(type, SID, 2, ip, port, reg,'joystick');
%% 3. RUN3
thermode_test(type, SID, 3, ip, port, reg,'joystick');
%% 4. RUN4
thermode_test(type, SID, 4, ip, port, reg,'joystick');
%% 5. RUN5
thermode_test(type, SID, 5, ip, port, reg,'joystick');
%% 6. RUN6
thermode_test(type, SID, 6, ip, port, reg,'joystick');













%% SEND DATA
SEMIC_Senddata(SID,'Calibration')