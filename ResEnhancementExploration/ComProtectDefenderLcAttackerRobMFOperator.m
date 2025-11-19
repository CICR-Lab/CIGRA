function [OptObj, BestDefendStrategy, BestAttackStrategy] = ComProtectDefenderLcAttackerRobMFOperator(System, DefenderParams, AttackerParams, OperatorParams, ComSet, varargin)
% ComProtectDefenderAttackerRobDCPFOperator
% --------------------------------------------------------
% Three-layer DAD model (Defender-Attacker-Defender) for MF systems
% Returns global optimal defense strategy and corresponding worst-case attack scenario
%
% Fixes: 
%   - Resolves iteration stalling issues
%   - Adds detailed debugging output
%   - Ensures all components in attack circle are included
%
% INPUT:
%   System          - System structure (with Node and Edge arrays)
%   DefenderParams  - Defense parameters
%   AttackerParams  - Attack parameters
%   OperatorParams  - System functionality evaluation parameters
%   ComSet          - Maximal component set
%   varargin        - Optional parameter: TerminalZone (for POP calculation)
%
% OUTPUT:
%   OptObj            - Optimal objective value (loss)
%   BestDefendStrategy- Structure {Node, Edge}
%   BestAttackStrategy- Structure {Node, Edge, Center}
%
% --------------------------------------------------------

% Process TerminalZone parameter
if nargin < 6
    TerminalZone = [];
else
    TerminalZone = varargin{1};
end

%% Step 1: Initialize bounds and variables
LowerObj = 0;           % Lower bound of objective value
UpperObj = inf;         % Upper bound of objective value
iteration = 0;          % Iteration counter
maxIterations = 10;     % Maximum number of iterations
tolerance = 1e-3;      % Convergence tolerance

% Initialize strategy storage
AttackStrategies = struct('Node', {}, 'Edge', {}, 'Center', {});
BestDefendStrategy = struct('Node', [], 'Edge', []);
BestAttackStrategy = struct('Node', [], 'Edge', [], 'Center', []);

% Get system size
N = length(System.Node);  % Number of nodes
E = length(System.Edge);  % Number of edges

%% Step 2: Iteratively solve DAD model
while UpperObj - LowerObj >= tolerance && iteration <= maxIterations
    iteration = iteration + 1;

    % === (1) Find worst-case attack under current defense strategy ===
    TempAttackerParams = AttackerParams;
    
    % Create invulnerable components vector (length N+E)
    TempAttackerParams.InvulnerableCom = zeros(1, N+E);
    
    % Set defended nodes
    if ~isempty(BestDefendStrategy.Node)
        % Ensure node indices are valid
        validNodes = BestDefendStrategy.Node(BestDefendStrategy.Node >= 1 & BestDefendStrategy.Node <= N);
        TempAttackerParams.InvulnerableCom(validNodes) = 1;
    end
    
    % Set defended edges
    if ~isempty(BestDefendStrategy.Edge)
        % Ensure edge indices are valid
        validEdges = BestDefendStrategy.Edge(BestDefendStrategy.Edge >= 1 & BestDefendStrategy.Edge <= E);
        TempAttackerParams.InvulnerableCom(N + validEdges) = 1;
    end
    
    TempAttackerParams.InvalidStrategy = [];

    % Call attacker function
    [AttackStrategyNow, AttackCenterNow, SysFunLoss, ~] = ...
        LcAttackerRobMFOperator(System, ComSet, OperatorParams, TempAttackerParams, TerminalZone);
    
    AttackStrategyNow = AttackStrategyNow.Unified;
    % Fix: Ensure attack strategy includes all components in attack circle
    % Find all components within attack circle
    allComponentsInCircle = [];
    for nc = 1:size(ComSet, 1)
        if norm(ComSet(nc, 1:3) - AttackCenterNow) < 1e-3
            temp = ComSet(nc, 4:end);
            temp = temp(temp >= 1 & temp <= N+E);
            allComponentsInCircle = temp; 
            break;
        end
    end
    
    % Add missing components to attack strategy
    if ~isempty(allComponentsInCircle)
        missingComponents = setdiff(allComponentsInCircle, AttackStrategyNow);
        if ~isempty(missingComponents)
            AttackStrategyNow = union(AttackStrategyNow, missingComponents);
        end
    end
    
    % Create attack strategy structure
    BestAttackStrategyNow = struct('Node', AttackStrategyNow(AttackStrategyNow <= N), ...
                                   'Edge', AttackStrategyNow(AttackStrategyNow > N) - N, ...
                                   'Center', AttackCenterNow);
    
    CurrentLoss = 1 - SysFunLoss(1);

    % Update lower bound if improved
    if LowerObj < CurrentLoss
        LowerObj = CurrentLoss;
    end

    % Check for convergence after attack optimization
    if UpperObj - LowerObj <= tolerance
        OptObj = LowerObj;
        BestAttackStrategy = BestAttackStrategyNow;
        break;
    end

    % === (1.1) Check for duplicate attack strategies ===
    InvalidStrategies = [];
    itag = 0;
    maxDuplicateChecks = 5; % Maximum attempts to find new strategy
    
    if ~isempty(AttackStrategies)
        while check_repeat_solution(AttackStrategies, BestAttackStrategyNow) == 1
            itag = itag + 1;
            
            % Break if max attempts reached
            if itag > maxDuplicateChecks
                break;
            end
            
            % Add to invalid strategies
            InvalidStrategies(itag).Node = BestAttackStrategyNow.Node;
            InvalidStrategies(itag).Edge = BestAttackStrategyNow.Edge;

            TempAttackerParams = AttackerParams;
            TempAttackerParams.InvulnerableCom = zeros(1, N+E);
            
            % Apply current defense strategy
            if ~isempty(BestDefendStrategy.Node)
                validNodes = BestDefendStrategy.Node(BestDefendStrategy.Node >= 1 & BestDefendStrategy.Node <= N);
                TempAttackerParams.InvulnerableCom(validNodes) = 1;
            end
            if ~isempty(BestDefendStrategy.Edge)
                validEdges = BestDefendStrategy.Edge(BestDefendStrategy.Edge >= 1 & BestDefendStrategy.Edge <= E);
                TempAttackerParams.InvulnerableCom(N + validEdges) = 1;
            end
            
            TempAttackerParams.InvalidStrategy = InvalidStrategies;

            % Find new attack strategy
            [AttackStrategyNow, AttackCenterNow, SysFunLoss, ~] = ...
                LcAttackerRobMFOperator(System, ComSet, OperatorParams, TempAttackerParams, TerminalZone);
            
            AttackStrategyNow = AttackStrategyNow.Unified;
            % Ensure all components in attack circle are included
            allComponentsInCircle = [];
            for nc = 1:size(ComSet, 1)
                if norm(ComSet(nc, 1:3) - AttackCenterNow) < 1e-3
                    temp = ComSet(nc, 4:end);
                    temp = temp(temp >= 1 & temp <= N+E);
                    allComponentsInCircle = temp;
                    break;
                end
            end
            
            % Add missing components
            if ~isempty(allComponentsInCircle)
                missingComponents = setdiff(allComponentsInCircle, AttackStrategyNow);
                if ~isempty(missingComponents)
                    AttackStrategyNow = union(AttackStrategyNow, missingComponents);
                end
            end
            
            % Create new attack strategy structure
            BestAttackStrategyNow = struct('Node', AttackStrategyNow(AttackStrategyNow <= N), ...
                                           'Edge', AttackStrategyNow(AttackStrategyNow > N) - N, ...
                                           'Center', AttackCenterNow);
        end
    end

    % Store new attack strategy
    AttackID = length(AttackStrategies);
    AttackStrategies(AttackID + 1).Node = BestAttackStrategyNow.Node;
    AttackStrategies(AttackID + 1).Edge = BestAttackStrategyNow.Edge;
    AttackStrategies(AttackID + 1).Center = BestAttackStrategyNow.Center;
    

    % === (2) Find optimal defense under current attack set ===
    DefenderParams.InvalidDefense = [];
    [UpperObj, BestDefendStrategy] = ComProtectDefenderGivenAttackStrategies(System, DefenderParams, AttackStrategies, OperatorParams);

    
    % Check for convergence after defense optimization
    if UpperObj - LowerObj <= tolerance
        OptObj = LowerObj;
        break;
    end
    
    % Break if max iterations reached
    if iteration == maxIterations
        OptObj = LowerObj;
        break;
    end
end

% Set optimal objective if not set in loop
if ~exist('OptObj', 'var')
    OptObj = LowerObj;
end

% Calculate final attack strategy
TempAttackerParams = AttackerParams;
TempAttackerParams.InvulnerableCom = zeros(1, N+E);

% Apply final defense strategy
if ~isempty(BestDefendStrategy.Node)
    validNodes = BestDefendStrategy.Node(BestDefendStrategy.Node >= 1 & BestDefendStrategy.Node <= N);
    TempAttackerParams.InvulnerableCom(validNodes) = 1;
end
if ~isempty(BestDefendStrategy.Edge)
    validEdges = BestDefendStrategy.Edge(BestDefendStrategy.Edge >= 1 & BestDefendStrategy.Edge <= E);
    TempAttackerParams.InvulnerableCom(N + validEdges) = 1;
end

% Find final attack strategy
[AttackStrategyNow, AttackCenterNow, ~, ~] = ...
    LcAttackerRobMFOperator(System, ComSet, OperatorParams, TempAttackerParams, TerminalZone);
AttackStrategyNow = AttackStrategyNow.Unified;

% Ensure all components in attack circle are included
allComponentsInCircle = [];
for nc = 1:size(ComSet, 1)
    if norm(ComSet(nc, 1:3) - AttackCenterNow) < 1e-3
        temp = ComSet(nc, 4:end);
        temp = temp(temp >= 1 & temp <= N+E);
        allComponentsInCircle = temp;
        break;
    end
end

% Add missing components
if ~isempty(allComponentsInCircle)
    missingComponents = setdiff(allComponentsInCircle, AttackStrategyNow);
    if ~isempty(missingComponents)
        AttackStrategyNow = union(AttackStrategyNow, missingComponents);
    end
end

% Create final attack strategy structure
BestAttackStrategy = struct('Node', AttackStrategyNow(AttackStrategyNow <= N), ...
                            'Edge', AttackStrategyNow(AttackStrategyNow > N) - N, ...
                            'Center', AttackCenterNow);

% Remove defended nodes from attack set
if ~isempty(BestDefendStrategy.Node) && ~isempty(BestAttackStrategy.Node)
    defendedNodes = BestDefendStrategy.Node;
    attackedNodes = BestAttackStrategy.Node;
    
    % Find intersection (defended nodes that are also attacked)
    conflictingNodes = intersect(defendedNodes, attackedNodes);
    if ~isempty(conflictingNodes)
        BestAttackStrategy.Node = setdiff(attackedNodes, defendedNodes);
    end
end

% Remove defended edges from attack set
if ~isempty(BestDefendStrategy.Edge) && ~isempty(BestAttackStrategy.Edge)
    defendedEdges = BestDefendStrategy.Edge;
    attackedEdges = BestAttackStrategy.Edge;
    
    % Find intersection (defended edges that are also attacked)
    conflictingEdges = intersect(defendedEdges, attackedEdges);
    if ~isempty(conflictingEdges)
        BestAttackStrategy.Edge = setdiff(attackedEdges, defendedEdges);
    end
end

% Display final results
fprintf('\n=== DAD Model Final Results ===\n');
fprintf('Iterations: %d\n', iteration);
fprintf('Optimal Objective: %.4f\n', OptObj);
fprintf('Defended nodes: %d\n', length(BestDefendStrategy.Node));
fprintf('Defended edges: %d\n', length(BestDefendStrategy.Edge));
fprintf('Attacked nodes: %d\n', length(BestAttackStrategy.Node));
fprintf('Attacked edges: %d\n', length(BestAttackStrategy.Edge));
fprintf('Attack center: [%.4f, %.4f, %.2f km]\n', BestAttackStrategy.Center(1), ...
                                          BestAttackStrategy.Center(2), ...
                                          BestAttackStrategy.Center(3));
end

%% Helper function: Check if attack strategy is repeated
function repeat_value = check_repeat_solution(AttackStrategies, BestAttackStrategy)
    repeat_value = 0;
    for s = 1:length(AttackStrategies)
        if isequal(sort(AttackStrategies(s).Node), sort(BestAttackStrategy.Node)) && ...
           isequal(sort(AttackStrategies(s).Edge), sort(BestAttackStrategy.Edge))
            repeat_value = 1;
            return;
        end
    end
end