function exp_scale(scale, joystick)

% explane scale (overall and continuous)
% THIS FUNCTION WILL RUN EVERY RUN.
% scale = 'predict'
% joystick = false or true, but it can be obtain own scripts.

global theWindow W H; % window property
global white red orange bgcolor; % color
global window_rect lb rb tb bb scale_H; % rating scale
global lb1 rb1 lb2 rb2;
global fontsize;


Screen(theWindow,'FillRect',bgcolor, window_rect);
velocity = cal_vel_joy('cont');
cir_center = [(lb1+rb1)/2 H*3/4+100];

switch scale
    case 'predict'
        msg = double('�� ô���� ȭ�鿡 ��Ÿ����\n "�̹� �ڱ��� �ִ� �󸶳� ���ñ��?"\n�� ���� ���� �� �ֽø� �ǰڽ��ϴ�.(Space)');
        % display
        Screen('TextSize', theWindow, fontsize);
        DrawFormattedText(theWindow, msg, 'center', 1/5*H, orange, [], [], [], 1.5);
        draw_scale('cont_predict_semicircular');
        Screen('Flip', theWindow);
        
        WaitSecs(1); %For preventing double-typed 'Space' key
        while (1)
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('space'))==1
                break
            elseif keyCode(KbName('q'))==1
                abort_experiment;
            end
        end
        %         Screen('TextSize', theWindow, fontsize);
        %         msg = double('�׷��� ���ݺ��� ������ �����ϰڽ��ϴ�.\n ô���� ȭ�鿡 �߸� �ٷ� ������ �ּ���.');
        %         DrawFormattedText(theWindow, msg, 'center', 'center', white, [], [], [], 1.5);
        %         Screen('Flip', theWindow);
        %        WaitSecs(1);
        rnd=randperm(2,1);
        WaitSecs(rnd);
        %
        %         cir_center = [(lb1+rb1)/2 H*3/4+100];
        %
        %         x=cir_center(1); y=cir_center(2);
        %         SetMouse(cir_center(1), cir_center(2)); % set mouse at the center
        %         start_while = GetSecs;
        %         while GetSecs-start_while < 11
        %             if joystick
        %                 [pos, ~] = mat_joy(0);
        %                 xAlpha=pos(1);
        %                 x=x+xAlpha*velocity;
        %                 yAlpha=pos(2);
        %                 y=y+yAlpha*velocity;
        %             else
        %                 [x,y,~]=GetMouse(theWindow);
        %             end
        %             % [x,y,~]=GetMouse(theWindow);
        %             radius = (rb1-lb1)/2; % radius
        %             theta = atan2(cir_center(2)-y,x-cir_center(1));
        %             % current euclidean distance
        %             curr_r = sqrt((x-cir_center(1))^2+ (y-cir_center(2))^2);
        %             % current angle (0 - 180 deg)
        %             curr_theta = rad2deg(-theta+pi);
        %             % For control a mouse cursor:
        %             % send to diameter of semi-circle
        %             if y > cir_center(2)%%bb
        %                 y = cir_center(2); %%bb;
        %                 SetMouse(x,y);
        %             end
        %             % send to arc of semi-circle
        %             if sqrt((x-cir_center(1))^2+ (y-cir_center(2))^2) > radius
        %                 x = radius*cos(theta)+cir_center(1);
        %                 y = cir_center(2)-radius*sin(theta);
        %                 SetMouse(x,y);
        %             end
        %             draw_scale('cont_predict_semicircular');
        %             Screen('DrawDots', theWindow, [x y], 15, orange, [0 0], 1);
        %             Screen('Flip', theWindow);
    case 'bar'
        x=cir_center(1)+randperm(100,1)-50; y=cir_center(2);
        radius = (rb2-lb2)/2;
        % Line
        theta = atan2(cir_center(2)-y,x-cir_center(1));
        x_arrow = radius*cos(theta)+cir_center(1);
        y_arrow = cir_center(2)-radius*sin(theta);
        xy = [cir_center(1) x_arrow; ...
            cir_center(2) y_arrow];
        
        
        msg = double('�� ô���� ȭ�鿡 ��Ÿ����\n "�̹� �ڱ��� �ִ� �󸶳� ���ñ��?"\n�� ���� 1) ���� 2) ��ȭ ���� �� �ֽø� �ǰڽ��ϴ�.(Space)');
        % display
        Screen('TextSize', theWindow, fontsize);
        DrawFormattedText(theWindow, msg, 'center', 1/5*H, orange, [], [], [], 1.5);
        
        draw_scale('overall_predict_semicircular');
        Screen(theWindow,'DrawLines', xy, 5, orange);
        Screen('Flip', theWindow);
        
        
        
               WaitSecs(1); %For preventing double-typed 'Space' key
        while (1)
            [~,~,keyCode] = KbCheck;
            if keyCode(KbName('space'))==1
                break
            elseif keyCode(KbName('q'))==1
                abort_experiment;
            end
        end
        %         Screen('TextSize', theWindow, fontsize);
        %         msg = double('�׷��� ���ݺ��� ������ �����ϰڽ��ϴ�.\n ô���� ȭ�鿡 �߸� �ٷ� ������ �ּ���.');
        %         DrawFormattedText(theWindow, msg, 'center', 'center', white, [], [], [], 1.5);
        %         Screen('Flip', theWindow);
        %        WaitSecs(1);
        rnd=randperm(2,1);
        WaitSecs(rnd);
end



end


