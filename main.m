clear

%% define graph parameters


% feel free to play around with the following parameters

graph_size = 100;
no_neighbours = 4;
rand_levels = 0:.5:1;

%% make graphs & plot them

[full_graph_regular, adjacency_graph_regular] = ...
    graph_generator(graph_size, no_neighbours, rand_levels(1));

[full_graph_sworld, adjacency_graph_sworld] = ...
    graph_generator(graph_size, no_neighbours, rand_levels(2));
   
[full_graph_random, adjacency_graph_random] = ...
    graph_generator(graph_size, no_neighbours, rand_levels(3));

% plot the graphs themselves
f1 = figure;

subplot(1,3,1)

plot(full_graph_regular,'Layout','circle','NodeColor','#F1BB7B', ...
    'EdgeColor','#F1BB7B');
title('Regular Graph')

subplot(1,3,2)
plot(full_graph_sworld,'Layout','circle','NodeColor', '#FD6467', ...
    'EdgeColor','#FD6467');
title('Small World Graph')

subplot(1,3,3)
plot(full_graph_random,'Layout','circle', 'NodeColor','#5B1A18', ...
    'EdgeColor','#5B1A18');
title('Random Graph')

% save the graph plots
saveas(gcf,'Graphs.png')


%% run Dijkstra & plot shortest distance descriptives

shortest_paths_regular = dijkstra(adjacency_graph_regular);
shortest_paths_sworld = dijkstra(adjacency_graph_sworld);
shortest_paths_random = dijkstra(adjacency_graph_random);

s_paths_regular_to_plot = shortest_paths_plottable(shortest_paths_regular);
s_paths_sworld_to_plot = shortest_paths_plottable(shortest_paths_sworld);
s_paths_random_to_plot = shortest_paths_plottable(shortest_paths_random);

% plot the shortest path distributions for each graph
f2 = figure;

subplot(1,3,1)

histogram(s_paths_regular_to_plot,'FaceColor','#F1BB7B');
title('Regular Graph Shortest Path Distribution')


subplot(1,3,2)
histogram(s_paths_sworld_to_plot,'FaceColor','#FD6467');
title('Small World Graph Shortest Path Distribution')


subplot(1,3,3)
histogram(s_paths_random_to_plot,'FaceColor','#5B1A18');
title('Random Graph Shortest Path Distribution')


% save the histograms
saveas(gcf,'Shortest_Path_Histogram.png')


% plot descriptive statistics
f3 = figure;

subplot(1,2,1)

boxchart([s_paths_regular_to_plot, s_paths_sworld_to_plot, ...
    s_paths_random_to_plot])

title('Boxplots')
xticklabels({'Regular Graph','Small World Graph','Random Graph'})


subplot(1,2,2)

means = [mean(s_paths_regular_to_plot), mean(s_paths_sworld_to_plot), ...
    mean(s_paths_random_to_plot)];
stds = [std(s_paths_regular_to_plot), std(s_paths_sworld_to_plot), ...
    std(s_paths_random_to_plot)];
errors = 2 * stds;

errorbar(means,errors,'-s','MarkerSize',10,...
    'MarkerEdgeColor','red','MarkerFaceColor','red');
xlim([0 4])
title('Mean \pm 2 Standard Deviations')
xticklabels({' ','Regular Graph','Small World Graph','Random Graph',' '})
 

sgtitle('Descriptive Statistics for Shortest Distances')

% save descriptive statistics
saveas(gcf,'Shortest_Path_Descriptives.png')


%% run Dijkstra on varying sizes of graphs & time performance

% this takes some time to run!
% feel free to adjust the candidate sizes for shorter run times
% but the long-ish runtimes are needed for meaningful comparisons

candidate_sizes = [500 1000 1500 2000];
candidate_neighbours = [10 20 50];
perf_rows = length(candidate_neighbours) * length(rand_levels);
perf_cols = length(candidate_sizes);

% below is the performance array where each column is
% a given candidate size and the rows are filled with
% combinations of candidate neighbours and randomness

performance = zeros(perf_rows,perf_cols);


% create a graph with ...
for c_size=1:length(candidate_sizes) % this size

    row_id = 1;

    for neighbour=1:length(candidate_neighbours) % this many neighbours

        row_id = row_id + 1;

        for randness=1:length(rand_levels) % this probability of randomness

            row_id = row_id + 1;

            % generate the graph
            [~,candidate_graph] = graph_generator(candidate_sizes(c_size),...
                candidate_neighbours(neighbour), rand_levels(randness));

            tic % start timing

            dijkstra(candidate_graph); % run the dijkstra algorithm

            performance(row_id, c_size) = toc; % stop timing and save time into array

        end
    end
end

% compute the mean and the standard deviation of performance
% for each graph size to plot

performance_means = mean(performance);
performance_errors = std(performance) *2;

f4 = figure;
errorbar(candidate_sizes,performance_means,performance_errors, ...
    '-s','MarkerSize',10,'MarkerEdgeColor','red','MarkerFaceColor','red');
title('Mean Runtime of Dijkstra Search for Different Graph Sizes')
xlabel('Number of Vertices')
ylabel('Time (s)')

saveas(gcf,'Dijkstra_Performance.png')