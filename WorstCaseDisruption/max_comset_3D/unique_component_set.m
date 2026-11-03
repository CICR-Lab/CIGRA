function ncom_set = unique_component_set(com_set, new_com_set)
% Inputs:
%   com_set - Existing component sets. Each row contains one attack center
%             (1-by-3) followed by component indices.
%   new_com_set - A candidate attack center and its component indices.
% Output:
%   ncom_set - Updated component sets with redundant combinations removed.

ntag = 0;

if ~isempty(com_set)
    ncom_set = com_set;
    [N, Nc] = size(com_set);        % N is the row count; Nc is the current width
    Nb = length(new_com_set);       % Nb is the candidate row width

    % Extract component indices after the first three center coordinates.
    coms = com_set(:, 4:Nc);        % Component indices in each existing row
    Bs = new_com_set(4:Nb);         % Component indices in the new set

    % Search the existing rows for equivalent component sets.
    temp = sum(ismember(coms, Bs), 2);           % Matches with Bs in each row
    num_non_zero = sum(coms > 0, 2);             % Valid components in each row
    mi = find(temp == num_non_zero);             % Rows whose valid components are all in Bs

    if ~isempty(mi)
        % Update the first duplicate row and remove any additional duplicates.
        ntag = 1;
        ncom_set(mi(1), 1:Nb) = new_com_set;
        if length(mi) > 1
            ncom_set(mi(2:end), :) = [];
        end
    else
        % Treat equal-size overlapping component sets as duplicates.
        if any(temp == Nb - 3)
            ntag = 1;
        end
    end

    % Add the candidate if it is a new component set.
    if ntag == 0
        if Nb <= Nc
            temp = zeros(1, Nc);
            temp(1:Nb) = new_com_set;
            ncom_set = [ncom_set; temp];
        else
            % Pad the existing matrix with zeros to reach Nb columns.
            ncom_set = [ncom_set, zeros(N, Nb - Nc); new_com_set];
        end
    end

else
    % Assign the candidate directly when the existing set is empty.
    ncom_set = new_com_set;
end
end
