"""Failure duration and magnitude analysis

This file includes both failure indices and calculates their correlation
using the Pearson correlation coefficient. The results are given both in graphical format
through the plot of correlation indices and in text format through the numerical output

CONTENTS:

- Evaluation of correlation between failure indices.
- Evaluation of correlation between elevation and failure indices.
- Evaluation of correlation between elevation and complex network metrics
- Normalization of CN, hydraulic and failures metrics
- Visual scatter plot of failure magnitude/duration correlation
- Visual scatter plot of failure duration and complex network metrics
- Visual scatter plot of failure magnitude and complex network metrics
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

networkList = [
    "c-town_true_network",  # 0
    "MOD",  # 1
    "KL",  # 2
    "ZJ",  # 3
    "Net3",  # 4
    "ky6",  # 5
    "Net1"  # 6
]

# Create a water network model 
inp_file = 'networks/' + networkList[1] + '.inp'
wn = wntr.network.WaterNetworkModel(inp_file)
print(inp_file)
water_link_diameter = wn.query_link_attribute('diameter')
water_link_length = wn.query_link_attribute('length')
invDia = 1 / water_link_diameter
node_elevation = wn.query_node_attribute('elevation')

# Creating graph from Water Network Model
G = wn.get_graph(node_weight=node_elevation, link_weight=water_link_length)
# G = wn.get_graph(wn)
sG = nx.Graph(G)
uG = G.to_undirected()  # undirected multigraph

# Adjust simulation time options for analyses
analysis_end_time = 7 * 24 * 3600
wn.options.time.duration = analysis_end_time
wn.options.time.report_timestep = 3600
wn.options.time.hydraulic_timestep = 3600  # 1hour

# Run a preliminary simulation to determine if junctions drop below the 
# pressure threshold during normal conditions

for name, junc in wn.junctions():
    junc.minimum_pressure = 0.0
    junc.nominal_pressure = 74.3
    Preq = junc.nominal_pressure

# sim = wntr.sim.WNTRSimulator(wn,mode = 'PDD')
sim = wntr.sim.WNTRSimulator(wn)
results = sim.run_sim()

# CNA metrics
betweenness_centrality = nx.betweenness_centrality(sG, weight='weight')
node_degree = dict(G.degree)
closeness_centrality = nx.closeness_centrality(G, distance='weight')
average_shortest_path_length = nx.average_shortest_path_length(uG, weight='weight')
print("APL:", round(average_shortest_path_length, 2))
print("Density:", round(nx.density(G), 3))
clustering_coeff = nx.clustering(sG, weight='weight')
eigenvector_centrality = nx.eigenvector_centrality_numpy(sG, weight='weight')
algebraic_connectivity = wntr.metrics.algebraic_connectivity(uG)
print("Algebraic connectivity:", round(algebraic_connectivity, 5))


# Evaluation of other CNA metrics, unrelated to this project purpose --> first stage exploration
bridges = wntr.metrics.bridges(G)
print("Density of bridges:", round(len(bridges) / G.number_of_edges(), 3))
central_point_dominance = wntr.metrics.central_point_dominance(G)
print("central point:", central_point_dominance)
print("Avg clustering coeff:", round(nx.average_clustering(sG, weight='weight'), 4))
mesh = (G.number_of_edges() - G.number_of_nodes()) / (2 * G.number_of_nodes() - 5)
print("Mesh:", mesh)


for name in wn.reservoir_name_list:  # Removing tank and reservoir from CNA metrics evaluatio
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

head = np.mean(results.node['head'])
pressure = results.node['pressure']

# Actual demand
actual_demand = results.node['demand'].loc[:, wn.junction_name_list]

# Actual demand sum for each node
actual_demand_tot = actual_demand.sum(axis=0)

# Expected demand
exp_demand = expected_demand(wn)
# print(exp_demand.to_string())

# Expected demand sum for each node
exp_demand_tot = exp_demand.sum(axis=0)

# Total demand all network
total_demand = exp_demand.values.sum()

# Exp tot - actual tot demand
sub = exp_demand_tot - actual_demand_tot
sub[sub < 0] = 0

# Failure duration index

sub_ex_ac = exp_demand.sub(actual_demand, axis=0)  # expected demand - actual demand
sub_ex_ac = sub_ex_ac.iloc[1:]  # exclude time = 0 in the evaluation, starting time is t = 1
tot_duration_fail = sub_ex_ac[sub_ex_ac > 0].count(axis=0)  # sum duration failure

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
                counter += 1
                failure += counter
                condition = False
                counter = 0
            elif column.values[i] <= 0:
                condition = True
                counter = 0
        tot_failure.append(failure)


failure_duration = tot_duration_fail.div(tot_failure).replace(np.NaN, 0)
# print("")
# print("Failure duration:")
# print(failure_duration)
# print("")
# print("Mean failure duration:")
# print(round(failure_duration.mean(),2))
# print("Median failure duration:")
# print(round(np.median(failure_duration), 2))

# Failure magnitude index
failure_magnitude = sub.div(total_demand)
failure_magnitude[failure_magnitude < 0] = 0
#print("")
# print("Mean failure magnitude:")
# print(round(failure_magnitude.mean(),6))
#print("Median failure magnitude")
#print(round(np.median(failure_magnitude), 6))

# Eigenvector centrality
eig = pd.Series(eigenvector_centrality)
normalized_eig = (eig - eig.min()) / (eig.max() - eig.min())

# Clustering coefficient
cluster_coeff = pd.Series(clustering_coeff)
normalized_cc = (cluster_coeff - cluster_coeff.min()) / (cluster_coeff.max() - cluster_coeff.min())

# Closeness centrality
close_cen = pd.Series(closeness_centrality)
normalized_close_cen = (close_cen - close_cen.min()) / (close_cen.max() - close_cen.min())

# Node degree
degree = pd.Series(node_degree)
normalized_degree = (degree - degree.min()) / (degree.max() - degree.min())

# BCE and failure
bce = pd.Series(betweenness_centrality)
normalized_bce = (bce - bce.min()) / (bce.max() - bce.min())

# Normalized failure durationd and magnitude
normalized_fd = (failure_duration - failure_duration.min()) / (failure_duration.max() - failure_duration.min())
normalized_fm = (failure_magnitude - failure_magnitude.min()) / (failure_magnitude.max() - failure_magnitude.min())

# Graphic plots

# Plot failure duration with color
# wntr.graphics.plot_network(wn, node_attribute=failure_duration)   # Darker colors correspond to higher failure values
# plt.title('Failure duration: '+ inp_file + ' ' + str(Preq))

# Plot failure magnitude color
# wntr.graphics.plot_network(wn, node_attribute=failure_magnitude)
# plt.title('Failure magnitude: '+ inp_file + ' ' + str(Preq))

# Correlation metrics

# Correlation Fail Mag/Fail Dur
print('Correlation failure duration/magnitude:', round(failure_duration.corr(failure_magnitude, method='pearson'), 3))

# ---- Correlation Elevation
print("---------------------------------------")
print("Correlation elevation/BCE:", round(node_elevation.corr(bce, method='pearson'), 3))
print("")
print("Correlation elevation/node_degree:", round(node_elevation.corr(degree, method='pearson'), 3))
print("")
print("Correlation elevation/Close_cen:", round(node_elevation.corr(close_cen, method='pearson'), 3))
print("")
print("Correlation elevation/Clustering Coefficient:", round(node_elevation.corr(cluster_coeff, method='pearson'), 3))
print("")
print("Correlation elevation/Eigenvector centrality:", round(node_elevation.corr(eig, method='pearson'), 3))
print("")
print("Correlation elevation/failure_duration:", round(node_elevation.corr(failure_duration, method='pearson'), 3))
print("")
print("Correlation elevation/failure_magnitude:", round(node_elevation.corr(failure_magnitude, method='pearson'), 3))


# Scatter plot magn/dur
# fig, ax = plt.subplots()
# ax.scatter(failure_duration, failure_magnitude, c = 'red', s = 10)
# plt.xlabel('Failure duration')
# plt.ylabel('Failure magnitude')
# plt.title('Correlation magnitude-duration: '+ inp_file + ' ' + str(Preq))


# Scatter plot fail duration
# fig, ax1 = plt.subplots()
# ax1.scatter(failure_duration, bce, c = 'red',alpha=0.4, s = 50)
# #ax.scatter(failure_duration, degree, c = 'yellow', marker = "x")
# ax1.scatter(failure_duration, close_cen, c = 'orange',alpha=0.5,s = 50)
# ax1.scatter(failure_duration, cluster_coeff, c = 'green', alpha=0.4, s = 50)
# ax1.scatter(failure_duration, eig, c = 'blue', marker = "^",alpha=0.4,s = 50)
# ax1.text(80, 0.6, 'BCE', color='red', fontsize=10)
# ax1.text(80, 0.57, 'CloseCen', color='orange', fontsize=10)
# ax1.text(80, 0.54, 'ClustCoef', color='green', fontsize=10)
# ax1.text(80, 0.51, 'EIG', color='blue', fontsize=10)
# plt.title('Scatter failure duration - CNA: '+ inp_file + ' ' + str(Preq))
# ===============================================================================


# Scatter plot normalized failure magnitude
fig, ax2 = plt.subplots()
ax2.scatter(normalized_fm, normalized_bce, c='red', alpha=0.4)
ax2.scatter(normalized_fm, normalized_degree, c='yellow', marker="x")
ax2.scatter(normalized_fm, normalized_close_cen, c='orange', alpha=0.4)
ax2.scatter(normalized_fm, normalized_cc, c='green', alpha=0.4)
ax2.scatter(normalized_fm, normalized_eig, c='blue', marker="^", alpha=0.4)
ax2.text(0.80, 0.97, 'BCE', color='red', fontsize=10)
ax2.text(0.80, 0.92, 'Degree', color='yellow', fontsize=10)
ax2.text(0.80, 0.87, 'CloseCen', color='orange', fontsize=10)
ax2.text(0.80, 0.82, 'ClustCoef', color='green', fontsize=10)
ax2.text(0.80, 0.77, 'EIG', color='blue', fontsize=10)
plt.title('Scatter failure magnitude - CNA: ' + inp_file + ' ' + str(Preq))

# Plot elevation
# wntr.graphics.plot_network(wn, node_attribute=node_elevation)
# plt.title('Elevation: '+ inp_file)


""" Todini index to evaluate water resilience

More info about Todini index:
<https://wntr.readthedocs.io/en/latest/apidoc/wntr.metrics.hydraulic.html?highlight=todini#wntr.metrics.hydraulic.todini_index> 

"""
results = sim.run_sim()   #Re-evaluation of simulation result for Todini index
head = results.node['head']
pressure = results.node['pressure']
demand = results.node['demand']
pump_flowrate = results.link['flowrate'].loc[:, wn.pump_name_list]
todini = wntr.metrics.todini_index(head, pressure, demand, pump_flowrate, wn, Pstar=Preq)
print("Todini:", np.mean(todini))

# plt.show()
