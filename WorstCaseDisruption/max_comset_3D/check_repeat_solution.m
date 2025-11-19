function repeat_value=check_repeat_solution(exist_solutions,current_solution)
repeat_value=0;
if ~isempty(exist_solutions)
   s=1;
   while s<=length(exist_solutions(1,:))
       if sum(exist_solutions(:,s)==current_solution)==length(current_solution)
          repeat_value=1;
          break;
       else
          s=s+1; 
       end
   end
end