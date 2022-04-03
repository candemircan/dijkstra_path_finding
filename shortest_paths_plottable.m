function shortest_paths = shortest_paths_plottable(shortest_paths)
% SHORTEST_PATHS rearranges the shortest path array to make it plottable

% >>> INPUT VARIABLES >>>
% NAME          TYPE            DESCRIPTION

% shortest_paths    matrix      matrix of shortest paths as 
%                               given by Dijkstra

% <<< OUTPUT VARIABLES <<<
% NAME          TYPE            DESCRIPTION

% shortest_paths    matrix      altered shortest path matrix for plotting
    
    shortest_paths = triu(shortest_paths); % only take the upper triangle
    shortest_paths = reshape(shortest_paths,1,[]); % flattent the matrix
    shortest_paths = nonzeros(shortest_paths); % remove the zeros

end