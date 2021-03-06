function rating_types = rating_prompt_repository

% rating_types = rating_prompt_repository(rating_types)
%
% Repository for rating types and prompts for the rating types


% ********** IMPORTANT NOTE **********
% YOU CAN ADD TYPES AND PROMPTS HERE. "cont_" AND "overall_" ARE IMPORTANT.
% * CRUCIAL: THE ORDER BETWEEN alltypes AND prompts SHOULD BE THE SAME.*
temp_rating_types = {'cont_int', '현재 느껴지는 통증의 세기를 지속적으로 보고해 주십시오.'; ...
    'cont_avoidance','앞으로 현재의 이 경험을 피하고 싶은 정도를 지속적으로 보고해주세요.'; ...
    'cont_avoidance_exp', '앞으로 현재의 이 경험을 피하고 싶은 정도를 지속적으로 보고해주세요.';...
    'overall_int', '자극이 얼마나 강했나요?';...
    'overall_avoidance', '이 경험을 피하고 싶은 정도를 보고해주세요.';...
    'overall_unpleasant', '본인이 경험한 불쾌함의 정도를 보고해주세요.'; ...
    'overall_aversive_ornot', '방금의 자극이 싫었나요?';...
    'overall_pain_ornot', '방금의 자극이 아팠나요?'; ...
    'overall_boredness', '방금 세션동안 얼마나 지겨우셨나요?'; ...
    'overall_alertness', '방금 세션동안 얼마나 졸리셨나요, 혹은 얼마나 정신이 또렸했나요?'; ...
    'overall_relaxed', '지금 얼마나 편안한 상태이신지요?';...
    'overall_attention', '방금 세션동안 과제에 얼마나 주의를 잘 기울이셨나요?'; ...
    'overall_resting_positive', '방금 세션동안 주로 했던 생각이 얼마나 긍정적이었나요?'; ...
    'overall_resting_negative', '방금 세션동안 주로 했던 생각이 얼마나 부정적이었나요?'; ...
    'overall_resting_myself', '방금 세션동안 주로 했던 생각이 나 자신에 대한 것이었나요?'; ...
    'overall_resting_others', '방금 세션동안 주로 했던 생각이 다른 사람들에 대한 것이었나요?'; ...
    'overall_resting_imagery', '방금 세션동안 주로 했던 생각이 생생한 이미지를 포함하고 있었나요?';...
    'overall_resting_present', '방금 세션동안 주로 했던 생각이 바로 지금, 여기에 대한 것이었나요?';...
    'overall_resting_past',  '방금 세션동안 주로 했던 생각이 과거에 대한 것이었나요?';...
    'overall_resting_future', '방금 세션동안 주로 했던 생각이 미래에 대한 것이었나요?';...
    'overall_resting_bitter_int', '방금 세션동안 경험했던 쓴 맛은 가장 심했을 때 얼마나 강했나요?';...
    'overall_resting_bitter_unp', '방금 세션동안 경험했던 쓴 맛은 가장 심했을 때 얼마나 불쾌했나요?';...
    'overall_resting_capsai_int', '방금 세션동안 경험했던 불타는 듯한 맛은 가장 심했을 때 얼마나 강했나요?';...
    'overall_resting_capsai_unp', '방금 세션동안 경험했던 불타는 듯한 맛은 가장 심했을 때 얼마나 불쾌했나요?';...
    'overall_thermal_int', '방금 세션동안 경험했던 열자극은 가장 심했을 때 얼마나 강했나요?';...
    'overall_thermal_unp', '방금 세션동안 경험했던 열자극은 가장 심했을 때 얼마나 불쾌했나요?';...
    'overall_pressure_int', '방금 세션동안 경험했던 압력자극은 가장 심했을 때 얼마나 강했나요?';...
    'overall_pressure_unp', '방금 세션동안 경험했던 압력자극은 가장 심했을 때 얼마나 불쾌했나요?';...
    'overall_negvis_int', '방금 세션동안 경험했던 시각자극은 가장 심했을 때 얼마나 강했나요?';...
    'overall_negvis_unp', '방금 세션동안 경험했던 시각자극은 가장 심했을 때 얼마나 불쾌했나요?';...
    'overall_negaud_int', '방금 세션동안 경험했던 소리자극은 가장 심했을 때 얼마나 강했나요?';...
    'overall_negaud_unp', '방금 세션동안 경험했던 소리자극은 가장 심했을 때 얼마나 불쾌했나요?';...
    'overall_posvis_int', '방금 세션동안 경험했던 시각자극은 가장 좋았을 때 얼마나 강했나요?';...
    'overall_posvis_ple', '방금 세션동안 경험했던 시각자극은 가장 좋았을 때 얼마나 기분이 좋았나요?';...
    'overall_mood', '방금 세션동안 기분은 어땠나요?';...
    'overall_comfortness', '지금 얼마나 편안한가요?';...
    'overall_avoidance_semicircular', '이 경험을 피하고 싶은 정도를 보고해주세요.';...
    'overall_motor', '목표 지점으로 움직인 후 클릭하세요.';...
    'overall_motor_semicircular', '목표 지점으로 움직인 후 클릭하세요.'};

rating_types.alltypes = temp_rating_types(:,1);
rating_types.prompts = temp_rating_types(:,2);

% korean to integer
for i = 1:numel(rating_types.prompts)
    rating_types.prompts{i} = double(rating_types.prompts{i});
end

end