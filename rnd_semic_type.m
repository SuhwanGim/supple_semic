function output = rnd_semic_type()
%%

%%
temp = ["dot","dot","dot","bar","bar","bar"];


while strcmp(temp,["dot","dot","dot","bar","bar","bar"]) == [1,1,1,1,1,1]
    rn = randperm(6,6); % 1 to 6)
    temp = temp(rn);
end

output = temp;


end


% 
% %%
% switch output(1)
%     case "bar"
%         disp('bar');
%     case "dot"
%         disp('dot');
%         
% end