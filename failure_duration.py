"""Failure duration analysis

Failure duration: Failure duration indicates how long a node is under failure and how long it takes
to recover, it determines if the state of a node is satisfactory (the pressure at that node is equal
or higher than the minimum requested one and water demand is fully supplied), or unsatisfactory,
e.g. the node pressure is under the minimum pressure threshold and the water demand cannot be fulfilled.
Failure duration is a node related measure and indicates the average duration time of failure occurred
to a specific node during simulation time.

CONTENTS:

- Calculation of failure duration index
- Calculation of complex network metrics
- Evaluation of correlation between head and failure magnitude indices.
- Evaluation of correlation between elevation and complex network metrics
- Visual scatter plot of complex network metrics
"""
import numpy as np
import networkx as nx
import wntr
from wntr.metrics.hydraulic import expected_demand, average_expected_demand
import matplotlib.pyplot as plt
from matplotlib import mlab
import pandas as pd
from wntr.metrics.topographic import algebraic_connectivity
from _operator import inv

# List of network (.inp file)
networkList = [
    "c-town_true_network",  #0
    "MOD",  #1
    "KL",  #2
    "ZJ",  #3
    "Net3",  #4
    "ky6",  #5
    "Net1"  #6
]

# Create a water network model
inp_file = 'networks/' + networkList[0] + '.inp'
wn = wntr.network.WaterNetworkModel(inp_file)
water_link_diameter = wn.query_link_attribute('diameter')
water_link_length = wn.query_link_attribute('length')
# invDia = 1/water_link_diameter
water_node_elevation = wn.query_node_attribute('elevation')

# Creating Graph from Water Network Model
# A graph can be tuned with or without weight

G = wn.get_graph(node_weight=water_node_elevation, link_weight=water_link_length)
# G = wn.get_graph(wn)
sG = nx.Graph(G)  # simple graph
uG = G.to_undirected()  # undirected multigraph

# Calculating Complex Network Metrics
betweenness_centrality = nx.betweenness_centrality(sG, weight='weight')
node_degree = dict(G.degree)
closeness_centrality = nx.closeness_centrality(G, distance='weight')
average_shortest_path_length = nx.average_shortest_path_length(uG, weight='weight')
clustering_coeff = nx.clustering(sG, weight='weight')
eigenvector_centrality = nx.eigenvector_centrality_numpy(sG, weight='weight')
algebraic_connectivity = wntr.metrics.algebraic_connectivity(uG)

for name in wn.reservoir_name_list:    # Removing tank and reservoir from CNA metrics evaluation
    del betweenness_centrality[name]
    del node_degree[name]
    del closeness_centrality[name]
    del clustering_coeff[name]
    del eigenvector_centrality[name]

for name in wn.tank_name_list:
    del betweenness_centrality[name]
    del node_degree[name]
    del closeness_centrality[name]
    del clustering_coeff[name]
    del eigenvector_centrality[name]

# Setting simulation time options for analyses
analysis_end_time = 24 * 7 * 3600
wn.options.time.duration = analysis_end_time
wn.options.time.report_timestep = 3600
wn.options.time.hydraulic_timestep = 3600  # 1 hour

# Run a preliminary simulation to determine if junctions drop below the
# pressure threshold during normal conditions

for name, junc in wn.junctions():
    junc.minimum_pressure = 0.0
    junc.nominal_pressure = 80
    Preq = junc.nominal_pressure

""" Setting simulation mode

- Selected PDD simulation

More info about pressure dependent demand (PDD) and demand-driven (DD) hydraulic simulation:
<https://wntr.readthedocs.io/en/latest/hydraulics.html> 

"""
sim = wntr.sim.WNTRSimulator(wn,mode = 'PDD')
# sim = wntr.sim.WNTRSimulator(wn) # DD mode
results = sim.run_sim()

# Water Network Metrics: Simulation results
head = np.mean(results.node['head'])
actual_demand = (results.node['demand'].loc[:, wn.junction_name_list])
exp_demand = expected_demand(wn)
pressure = results.node['pressure']

# print(pressure.max(axis = 1))
# print("Mediana 10 percentile:", round(med10,1))
# print("Mediana 50 percentile:", round(med50,1))
# print("Mediana 90 percentile:", round(med90,1))


# Plotting Empirical Cumulative Distribution Function

# def ecdf(data):
#     """ Compute ECDF """
#     x = np.sort(data)
#     n = x.size
#     y = np.arange(1, n+1) / n
#     return(x,y)
#     
# x,y = ecdf(pressure.values.tolist())
# plt.grid(True)
# plt.title("ECDF")
# plt.xlabel("Pressure")
# med10 = np.percentile(pressure, 10)
# med50 = np.median(pressure)
# med90 = np.percentile(pressure, 90)
# ind = [0, med10 ,med50 ,med90]
# plt.xticks(ind, ('0','10','50','90'))
# plt.scatter(x=x,y=y, s = 10, c = 'red', alpha = 0.6)
# plt.show()


# Evaluating Failure Duration

sub_ex_ac = exp_demand.sub(actual_demand, axis=0)  # expected demand - actual demand
sub_ex_ac = sub_ex_ac.iloc[1:]  # exclude time = 0 in the evaluation, starting time is t = 1
tot_duration_fail = sub_ex_ac[sub_ex_ac > 0].count(axis=0)  # sum duration failure

# print("Sub expected - actual:")
# print(sub_ex_ac.to_string())
# print(sum(sum_fail.axes[0] > 0))
# print(sub_ex_ac[sub_ex_ac > 0].index)


tot_failure = []

for index, column in sub_ex_ac.iteritems():
    condition = True;
    counter = 0
    failure = 0
    if (column == 0).all():
        tot_failure.append(0)
    elif (column > 0).all():
        tot_failure.append(1)
    else:
        for i in range(len(column.values)):
            counter = 0
            if column.values[i] > 0 and condition == True:
                counter = counter + 1
                failure = failure + counter
                condition = False
                counter = 0
            elif column.values[i] <= 0:
                condition = True
                counter = 0
        tot_failure.append(failure)

print("")
# print("Number of total failure:",*tot_failure)


# Failure duration index

failure_duration = tot_duration_fail.div(tot_failure).replace(np.NaN, 0) # If the value is not valid of negative, replace it with 0
# print("")
# print("Failure duration:")
# print(failure_duration)
# print("")
# print("Mean failure duration:" + round(np.median(failure_duration),2))

# Plot elevation
wntr.graphics.plot_network(wn, node_attribute=water_node_elevation)
plt.title("Elevation")

# Plot water network based on the severity of the failure duration
wntr.graphics.plot_network(wn, node_attribute=failure_duration)  # Darker colors correspond to higher failure values
plt.title("Failure duration")

# CNA metrics plot and their correlation with failure magnitude

# Betweenness centrality
bce = pd.Series(betweenness_centrality)
wntr.graphics.plot_network(wn, node_attribute=betweenness_centrality)
# plt.title("BCE")
print("")
print("Correlation Failure duration/BCE:", round(failure_duration.corr(bce, method='pearson'), 3))

# ===============================================================================
# fig, ax = plt.subplots()
# ax.scatter(failure_duration, bce, c = 'red')
# ax.plot(np.unique(failure_duration), np.poly1d(np.polyfit(failure_duration, bce, 1))(np.unique(failure_duration)), color='yellow')
# plt.xlabel('Failure duration')
# plt.ylabel('BCE')
# ===============================================================================


# Node degree
degree = pd.Series(node_degree)
# wntr.graphics.plot_network(wn, node_attribute=degree)
# plt.title("Node degree")
print("")
print("Correlation Failure duration/node degree:", round(failure_duration.corr(degree, method='pearson'), 3))

# ===============================================================================
# fig, ax = plt.subplots()
# ax.scatter(failure_duration, degree, c = 'red')
# ax.plot(np.unique(failure_duration), np.poly1d(np.polyfit(failure_duration, degree, 1))(np.unique(failure_duration)), color='yellow')
# plt.xlabel('Failure duration')
# plt.ylabel('degree')
# ===============================================================================

# Closeness centrality
close_cen = pd.Series(closeness_centrality)
# wntr.graphics.plot_network(wn, node_attribute=close_cen)
# plt.title("Closeness centrality")
print("")
print("Correlation Failure duration/Closeness Centrality:",
      round(failure_duration.corr(close_cen, method='pearson'), 3))

# ===============================================================================
# fig, ax = plt.subplots()
# ax.scatter(failure_duration, close_cen, c = 'red')
# ax.plot(np.unique(failure_duration), np.poly1d(np.polyfit(failure_duration, close_cen, 1))(np.unique(failure_duration)), color='yellow')
# plt.xlabel('Failure duration')
# plt.ylabel('Closeness')
# ===============================================================================


# Clustering coefficient
cluster_coeff = pd.Series(clustering_coeff)
# wntr.graphics.plot_network(wn, node_attribute=cluster_coeff)
# plt.title("Clustering coefficient")
print("")
print("Correlation Failure duration/Clustering Coefficient:",
      round(failure_duration.corr(cluster_coeff, method='pearson'), 3))

# ===============================================================================
# fig, ax = plt.subplots()
# ax.scatter(failure_duration, cluster_coeff, c = 'red')
# ax.plot(np.unique(failure_duration), np.poly1d(np.polyfit(failure_duration, cluster_coeff, 1))(np.unique(failure_duration)), color='yellow')
# plt.xlabel('Failure duration')
# plt.ylabel('Clustering coefficient')
# ===============================================================================

# Eigenvector centrality
eig = pd.Series(eigenvector_centrality)
# wntr.graphics.plot_network(wn, node_attribute=cluster_coeff)
# plt.title("EIG centrality")
print("")
print("Correlation Failure duration/Eigenvector centrality:",
      round(failure_duration.corr(eig, method='pearson'), 3))

# ===============================================================================
# fig, ax = plt.subplots()
# ax.scatter(failure_duration, eig, c = 'red')
# ax.plot(np.unique(failure_duration), np.poly1d(np.polyfit(failure_duration, eig, 1))(np.unique(failure_duration)), color='yellow')
# plt.xlabel('Failure duration')
# plt.ylabel('EIG')
# ===============================================================================


# Correlation Head with failure duration
print("")
print("Correlation head/failure_duration:", round(head.corr(failure_duration, method='pearson'), 3))

# Correlation Elevation with CNA metrics and failure duration
print("---------------------------------------")
print("Correlation elevation/BCE:", round(water_node_elevation.corr(bce, method='pearson'), 3))
print("")
print("Correlation elevation/node_degree:", round(water_node_elevation.corr(degree, method='pearson'), 3))
print("")
print("Correlation elevation/Close_cen:", round(water_node_elevation.corr(close_cen, method='pearson'), 3))
print("")
print("Correlation elevation/Clustering Coefficient:", round(water_node_elevation.corr(cluster_coeff, method='pearson'), 3))
print("")
print("Correlation elevation/Eigenvector centrality:", round(water_node_elevation.corr(eig, method='pearson'), 3))
print("")
print("Correlation elevation/failure_duration:", round(water_node_elevation.corr(failure_duration, method='pearson'), 3))

# ---- Scatter plot ----#
fig, ax = plt.subplots()
ax.scatter(failure_duration, bce, c='red')
# ax.scatter(failure_duration, degree, c = 'yellow', marker = "x")
ax.scatter(failure_duration, close_cen, c='green', marker="2")
ax.scatter(failure_duration, cluster_coeff, c='blue', marker="+")
ax.scatter(failure_duration, eig, c='pink', marker="^")

# plt.show()
