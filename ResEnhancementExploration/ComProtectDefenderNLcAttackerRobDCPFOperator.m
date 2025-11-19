function [BestDefendStrategy,BestAttackStrategy,OptObj]=ComProtectDefenderNLcAttackerRobDCPFOperator(System,DefenderParams,AttackerParams,OperatorParams)
% rob_defender_attacker_operator_COMP_NPCA_MF solves the max-min-max objective_function problem, returns the optimal objective, best defense strategy
% and the worst-case attack strategy under the best defense
%
% input:
%     System              struct with fields .Node and .Edge (power, gas or water)
%       System.Node(n).ID
%       System.Node(n).RealDemand, .TargetDemand, .RealGeneration, .MaxGeneration
%       System.Edge(e).ID
%       System.Edge(e).FromNodeID, .ToNodeID, .Length, .RealFlow, .Capacity
%
%     DefenderParams        struct
%       .Budget             scalar total defend budget
%       .NodeDedendCost     N×1 node costs
%       .EdgeDefendCost     E×1 edge costs
%       .InvulNode          IDs of invulnerable nodes
%       .InvulEdge          IDs of invulnerable edges
%       .InvalidStrategy(k)  struct with .Node and .Edge lists
%
%     AttackerParams        struct
%       .Budget             scalar total attack budget
%       .NodeAttackCost     N×1 node costs
%       .EdgeAttackCost     E×1 edge costs
%       .InvulNode          IDs of invulnerable nodes
%       .InvulEdge          IDs of invulnerable edges
%       .InvalidStrategy(k)  struct with .Node and .Edge lists
%
%     OperatorParams      struct
%       .NodeWeight         N×1 node‐weight (defaults to ones(N,1))
%
%   Outputs
%     OptObj              scalar, the minimal attack consequence
%     BestDefendStrategy  struct
%       .Node             IDs of defended nodes
%       .Edge             IDs of defended edges
%    
%     BestAttackStrategy  struct
%       .Node             IDs of attacked nodes
%       .Edge             IDs of attacked edges

%% set model parameters and initialize solutions
LowerObj=0;
UpperObj=inf;
xhcs=0;

%% define those identified attack strategies and also the best defend strategy
AttackStrategies=[];
BestDefendStrategy=[];

while UpperObj-LowerObj>=0.0001

    xhcs=xhcs+1;

    [xhcs LowerObj UpperObj]

    % solve the best attack strategy under given a defense strategy
    TempAttackerParams=AttackerParams;
    if ~isempty(BestDefendStrategy)
        TempAttackerParams.InvulNode=BestDefendStrategy.Node;
        TempAttackerParams.InvulEdge=BestDefendStrategy.Edge;
    end
    TempAttackerParams.InvalidStrategy=[];

    [BestAttackStrategy,~,temp_lower_obj]=NLcAttackerRobDCPFOperator(System,TempAttackerParams,OperatorParams);
    if LowerObj<temp_lower_obj
        LowerObj=temp_lower_obj;
    end

    if UpperObj-LowerObj<=0.0001
        OptObj=LowerObj;
        break;
    end

    InvalidStrategies=[];itag=0;
    % check whether the solved best attack strategy is repeated with previous solutions,if repeated, need to update the best attack strategy
    if ~isempty(AttackStrategies)
        while check_repeat_solution(AttackStrategies,BestAttackStrategy)==1
            itag=itag+1;
            InvalidStrategies(itag).Node=BestAttackStrategy.Node;
            InvalidStrategies(itag).Edge=BestAttackStrategy.Edge;

            TempAttackerParams=AttackerParams;
            if ~isempty(BestDefendStrategy)
                TempAttackerParams.InvulNode=BestDefendStrategy.Node;
                TempAttackerParams.InvulEdge=BestDefendStrategy.Edge;
            end
            TempAttackerParams.InvalidStrategy=InvalidStrategies;

            [BestAttackStrategy,~]=NLcAttackerRobDCPFOperator(System,TempAttackerParams,OperatorParams);
        end
    end
    AttackID=length(AttackStrategies);
    AttackStrategies(AttackID+1).Node=BestAttackStrategy.Node;
    AttackStrategies(AttackID+1).Edge=BestAttackStrategy.Edge;

    % solve the best defense under given attack strategies
    DefenderParams.InvalidDefense=[];
    [UpperObj,BestDefendStrategy]=ComProtectDefenderGivenAttackStrategies(System,DefenderParams,AttackStrategies,OperatorParams);

    if UpperObj-LowerObj<=0.0001
        OptObj=LowerObj;

        TempAttackerParams=AttackerParams;
        if ~isempty(BestDefendStrategy)
            TempAttackerParams.InvulNode=BestDefendStrategy.Node;
            TempAttackerParams.InvulEdge=BestDefendStrategy.Edge;
        end
        TempAttackerParams.InvalidStrategy=InvalidStrategies;

        [BestAttackStrategy,~]=NLcAttackerRobDCPFOperator(System,TempAttackerParams,OperatorParams);
        break;
    end
end

    function repeat_value=check_repeat_solution(AttackStrategies,BestAttackStrategy)
        repeat_value=0;
        if ~isempty(AttackStrategies)
            s=1;
            while s<=length(AttackStrategies)
                if length(intersect(AttackStrategies(s).Node,BestAttackStrategy.Node))==length(BestAttackStrategy.Node) && length(intersect(AttackStrategies(s).Edge,BestAttackStrategy.Edge))==length(BestAttackStrategy.Edge)
                    repeat_value=1;
                    break;
                else
                    s=s+1;
                end
            end
        end
    end
end