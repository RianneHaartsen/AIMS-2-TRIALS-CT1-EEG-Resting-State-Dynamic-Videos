function [ indices ] = findcontig2( search, sought )
    %UNTITLED2 Summary of this function goes here
    %   Detailed explanation goes here

    indices=[];
    
    if ~isvector(search)
        error('Search variable must be a vector')
    end
    
    % ensure search is a col vector
    if size(search, 1) > size(search, 2)
        search = search';
    end
    
    % if sought is not specified, assume true
    if ~exist('sought', 'var') || isempty(sought)
        sought = true;
    end

    found = search==sought;
    if ~any(found)
        % sought value not found in search array
        return
    elseif all(found)
        % search array was full of sought values
        indices = [1, length(search), length(search)];
        return
    end
    
    lag1 = [found, false];
    lag2 = [false, found];

    % find onset and offsets
    onset = lag1 & ~lag2;
    offset = ~lag1 & lag2;
    
    % find indices of on/offsets
    idxOnset = find(onset);
    idxOffset = find(offset) - 1;

    % calculate duration, make table
    indices = [idxOnset', idxOffset', idxOffset' - idxOnset' + 1];

end

