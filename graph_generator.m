function [my_graph,adjacency_mat]= graph_generator(vertices, neighbours, randomness)
% GRAPH_GENERATOR generates graphs with the given parameters
%   The procedure is implemented in the way as described in https://www.nature.com/articles/30918

% >>> INPUT VARIABLES >>>
% NAME          TYPE            DESCRIPTION

% vertices      scalar          total number of vertices of the graph
% neighbours    scalar          initial number of nearest neighbours
%                               a vertex connects with
% randomness    scalar          the probability that a given edge
%                               will reconnect with a different vertex

% <<< OUTPUT VARIABLES <<<
% NAME          TYPE            DESCRIPTION

% my_graph      graph           graph defined by the final adjacancy_mat
% adjacancy_mat 2-D matrix      vertices x vertices adjacency matrix

    
% ------- checks of input arguments ---------------------------------

    if randomness < 0 || randomness > 1
        error("Randomness must be between 0 & 1.")
    end

    if neighbours >= vertices
        error("Number of neighbours can at most be edges - 1.")
    end

    if ~(mod(neighbours,1) == 0 && mod(vertices,1) == 0)
        error("Neighbours and edges must be integers")
    end

% ------- end of input checks --------------------------------------
    
    
    
    % because we are working with undirected & unweighted graphs here,
    % our adjacancy matrix (A) is a matrix with binary entries where
    % the entry A(i,j) defines whether the vertices i & j are
    % connected (==1) or not (==0)

    adjacency_mat = zeros(vertices); % initialise the adjacancy matrix

    % create an adjacancy row which defines how a given vertex is
    % connected to other vertices. In this case, the vertex's n/2
    % nearest neighbours on each side should be 1, and all other entries
    % should be zero
    
    adjacency_row = zeros(1,vertices); 
    neighbour_set_1 = floor(neighbours/2);
    adjacency_row(1:neighbour_set_1) = 1;
    neighbour_set_2 = neighbours - neighbour_set_1 -1 ;
    start_point_2 = neighbour_set_1+2;
    adjacency_row(start_point_2: start_point_2 + neighbour_set_2) = 1;
    adjacency_row = circshift(adjacency_row,-(start_point_2-2));

    % conveniently define adjacancy matrix by sliding the adjacency row
    % forward as we progress through the rows of the matrix

    for vertex=1:vertices
        adjacency_mat(vertex,:) = adjacency_row;
        adjacency_row = circshift(adjacency_row,1);
    end

    % unweighted & undirected graphs are symmetric
    % make the adjacancy matrix symmetric
    
    adjacency_mat = make_symmetric(adjacency_mat,"upper");

    % not all combinations of vertices and nearest neighbours
    % can be made (for example a 3 vertex graph with 1 nearest
    % neighbour cannot be symmetric and must have one vertex
    % with no edges). Check if the given combination is possible
    % otherwise carry on with the symmetric matrix created above
    % as it is a close approximation. Warn the user about it!
    
    if any(sum(adjacency_mat,2) ~= neighbours * ones(vertices,1))
        warning(['A symmetric graph of size %d with %d '...
        ,'nearest neighbour connections is not possible', ...
        '\n This function will return a close approximation'],...
        round(vertices),round(neighbours))
    end
    
 
    % if the randomness parameter is > 0
    if randomness
        for row=1:(vertices-1) % last vertex is defined by the other vertices' edges
            for col=1:(vertices-1)
                if reassign() % reconnect the edge if multiple conditions are met
                    adjacency_mat(row,col) = 0; % remove the existing edge
                    new_neighbour = randi([row+1 vertices]); % choose the new edge randomly
                    adjacency_mat(row, new_neighbour) = 1; % add the new edge

                end
            end
        end
    end


    adjacency_mat = make_symmetric(adjacency_mat,"upper"); %make the final graph symmetric
    my_graph = graph(adjacency_mat); %turn the adjacancy matrix into a graph for plotting

    function reassign_this = reassign()
    % REASSIGN_THIS checks for the multiple conditions needed to reassign an edge
    % note that this is a nested function only used by its parent function

    % <<< OUTPUT VARIABLES <<<
    % NAME          TYPE            DESCRIPTION

    % reassign_this logical         whether reassignment of edge is allowed or not

        reassign_this = false; % initialise the reassignment allowment to false
        if row > col % if we are in the upper half of the matrix
            if sum(rand >= cumsum([1-randomness, randomness])) % reassignment with probability randomness
                if adjacency_mat(row,col) == 1 % if there already exists an edge between the vertices
    
                    reassign_this = true;
                
                end
            end
        end
    end

end

function symmetric_matrix = make_symmetric(some_matrix,to_copy)
% SYMMETRIC_MATRIX makes the given matrix symmetric
% note that this function can only be used within this .m file

% >>> INPUT VARIABLES >>>
% NAME          TYPE            DESCRIPTION

% some_matrix   2-D matrix       input array
% to_copy       char            'upper' copies the upper triangle to lower
%                               'lower' copies the lower triangle to upper

% <<< OUTPUT VARIABLES <<<
% NAME          TYPE            DESCRIPTION

% symmetric_matrix  2-D matrix  final symmetric array

    
% ------- checks of input arguments ---------------------------------
    
    if (to_copy ~="lower" && to_copy ~="upper")
        error('to_copy must be either lower or upper')
    end
    
    if ~ismatrix(some_matrix)
        error('The first input must be a matrix')
    end
% ------- checks of input arguments ---------------------------------    
    
    
    switch to_copy
        case 'upper'
            symmetric_matrix = triu(some_matrix) + triu(some_matrix,1)';
        case 'lower'
            symmetric_matrix = tril(some_matrix) + tril(some_matrix,-1)';
    end
end