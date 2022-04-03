function shortest_paths = dijkstra(my_graph)
% DIJKSTRA computes the shortest path between all edges of the graph

% >>> INPUT VARIABLES >>>
% NAME          TYPE            DESCRIPTION

% my_graph      matrix          adjacancy matrix defining a graph

% <<< OUTPUT VARIABLES <<<
% NAME          TYPE            DESCRIPTION

% shortest_paths matrix         a matrix M where entry M(i,j)
%                               defines the shortest path between
%                               the points i & j


    [~,nodes] = size(my_graph); % get graph size
    shortest_paths = inf(nodes); % initialise all shortest paths to max values
    for start_point=1:nodes % from the start point...
        shortest_path_source = inf(1,nodes); % shortest path from start_point to vertices
        shortest_nodes = zeros(1,nodes); % initialise the shortest visited vertices from start_point to false
        shortest_path_source(start_point) = 0; % set distance to self to 0
        for step=1:nodes
            min_dist = inf; % initial minimum distance is infinity
            for vertex=1:nodes % check if any of the unvisited vertices are of closer distance
                if (shortest_path_source(vertex) < min_dist && ...
                        ~shortest_nodes(vertex))
                    min_dist = shortest_path_source(vertex); % record the closer distance
                    min_dist_idx = vertex; % record the closer vertex
                end
            end
    
            shortest_nodes(min_dist_idx) = true; % mark the closest vertex
    
            for end_point=1:nodes % check for each target vertex...
                if (my_graph(min_dist_idx,end_point) > 0 && ... % that it is different from the start vertex
                        ~shortest_nodes(end_point) && ... % that it is not the closest vertex
                        shortest_path_source(end_point) > ... % that the currently recorded distance ...
                        shortest_path_source(min_dist_idx) + ... % is bigger than the distance to be computed ...
                        my_graph(min_dist_idx,end_point)) % by adding shortest distance to the edge information
                    
                    shortest_path_source(end_point) = shortest_path_source(min_dist_idx) + ...
                        my_graph(min_dist_idx,end_point); % add edge to the shortest distance to compute the distance to end_point
    
                end
            end
    
        end
        shortest_paths(start_point,:) = shortest_path_source; % record shortest_path_source as a row in the shortest_paths matrix
    end



end