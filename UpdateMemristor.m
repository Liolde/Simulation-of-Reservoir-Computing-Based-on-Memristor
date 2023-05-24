function [update_G] = UpdateMemristor(G_val,bool,rise,decline)

if(bool)
    [~,index] = min(abs(rise-G_val));
    update_G = rise(index+1000);
else
    [~,index] = min(abs(decline-G_val));
    update_G = decline(index+1000);
end

