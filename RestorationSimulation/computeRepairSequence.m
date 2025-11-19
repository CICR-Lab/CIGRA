function res_step = computeRepairSequence(res_seq, RR)
% Convert priority list (res_seq) into a chronologic schedule with RR teams.
%
% INPUT:
%    res_seq  - [DamageType, ComponentID, RepairTime, SysType]
%    RR       -  number of repair teams
% OUTPUT:
%    res_step - [DamageType, ComponentID, SysType, FinishTime, TeamID]

Ndm = size(res_seq,1);
if Ndm==0
    res_step=[];
    return
end
res_step = zeros(Ndm,5);

if RR == 1
    finish_time = cumsum(res_seq(:,3));  % sum of RepairTime
    res_step(:,1) = res_seq(:,1);        % DamageType
    res_step(:,2) = res_seq(:,2);        % ComponentID
    res_step(:,3) = res_seq(:,4);        % SysType
    res_step(:,4) = finish_time;         % FinishTime
    res_step(:,5) = 1;                   % TeamID
else
    team_available = zeros(RR,1); % earliest available time for each team
    for i = 1:Ndm
        duration = res_seq(i,3);
        [avail_time, team_id] = min(team_available);
        finish_time = avail_time + duration;
        team_available(team_id) = finish_time;
        res_step(i,:) = [res_seq(i,1), res_seq(i,2), res_seq(i,4), finish_time, team_id];
    end
    % sort chronologically by finish time
    res_step = sortrows(res_step, 4);
end
end