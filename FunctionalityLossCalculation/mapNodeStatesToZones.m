function ZoneState = mapNodeStatesToZones(CIS, nodeState, TerminalZone, Type)
% mapNodeStatesToZones
%   Map node-level functionality states at a single time slice to
%   zone-level states based on CIS.Node(n).ServiceZone(:,1).
%
% INPUTS:
%   CIS        每 CIS struct with Node(i).ServiceZone containing zone indices.
%   nodeState  每 N℅1 vector of node-level states at one time slice
%                (e.g., pre-disaster or post-disaster).
%   TerminalZone 每 structure array of zones; only used here to define the
%                total number of zones (numel(TerminalZone)).
%   Type         - Node state type: 'Flow' or 'Topology'
% OUTPUT:
%   ZoneState  每 G℅1 vector of zone-level states

numZones  = numel(TerminalZone);
ZoneState = nan(numZones, 1);
N = numel(CIS.Node);

switch Type
    case 'Flow'
        Denom = [CIS.Node.TargetDemand]';
        F_node = max(0, min(1, nodeState./Denom));
        
        for n = 1:N
            if ~isempty(CIS.Node(n).ServiceZone)
                g = CIS.Node(n).ServiceZone(:,1);
                ZoneState(g) = F_node(n);
            end
        end
    case 'Topology'
        for n = 1:N
            if ~isempty(CIS.Node(n).ServiceZone)
                g = CIS.Node(n).ServiceZone(:,1);
                ZoneState(g) = nodeState(n);
            end
        end
end
end