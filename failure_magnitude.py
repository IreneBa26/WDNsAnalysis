"""Failure magnitude analysis

Failure magnitude: evaluates how bad the water shortage is at each node, it is expressed
as the percentage of unmet water demand. Results are aggregated over all time-steps.

CONTENTS:

- Calculation of failure magnitude index
- Calculation of complex network metrics
- Evaluation of correlation between head and failure magnitude indices.
- Visual scatter plot of complex network metrics
"""

import networkx as nx
import wntr
import matplotlib.pyplot as plt
import winsound
import numpy as np
import pandas as pd
from wntr.network.elements import Junction, Demands
from numpy import short
from wntr.morph.link import break_pipe
import wntr.network.controls as controls
from wntr.metrics.hydraulic import expected_demand, average_expected_demand
from wntr.epanet.util import FlowUnits, HydParam
from wntr.epanet.util import *
from networkx.algorithms.bipartite.centrality import betweenness_centrality

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

inp_file = 'networks/' + networkList[6] + '.inp'  # Select network from list
print(inp_file)
wn = wntr.network.WaterNetworkModel(inp_file)
water_link_diameter = wn.query_link_attribute('diameter')
invDia = 1 / water_link_diameter
water_link_length = wn.query_link_attribute('length')
node_elevation = wn.query_node_attribute('elevation')

# Creating graph from Water Network Model
G = wn.get_graph(node_weight=node_elevation, link_weight=water_link_length)
# G = wn.get_graph(wn)
sG = nx.Graph(G)
uG = G.to_undirected()  # undirected multigraph

# Complex Network Analysis metrics
betweenness_centrality = nx.betweenness_centrality(sG, weight='weight')
node_degree = dict(G.degree)
closeness_centrality = nx.closeness_centrality(G, distance='weight')
clustering_coeff = nx.clustering(sG, weight='weight')
eigenvector_centrality = nx.eigenvector_centrality_numpy(sG, weight='weight')

for name in wn.reservoir_name_list:   # Removing tank and reservoir from CNA metrics evaluation
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

# Setting simulation time parameters
wn.options.time.duration = 7 * 24 * 3600  # week 7*24*3600
wn.options.time.report_timestep = 3600
wn.options.time.hydraulic_timestep = 3600  # 1hour

# Setting minimum e requested pressure
for name, junc in wn.junctions():
    # print(junc.demand_timeseries_list[0])
    # junc.demand_timeseries_list[0].base_value = junc.demand_timeseries_list[0].base_value*2
    junc.minimum_pressure = 0.0
    junc.nominal_pressure = 80

""" Setting simulation mode

- Selected PDD simulation

More info about pressure dependent demand (PDD) and demand-driven (DD) hydraulic simulation:
<https://wntr.readthedocs.io/en/latest/hydraulics.html> 

"""

sim = wntr.sim.WNTRSimulator(wn, mode='PDD')
results = sim.run_sim()
head = results.node['head']
pressure = results.node['pressure']
demand = results.node['demand']
pump_flowrate = results.link['flowrate'].loc[:, wn.pump_name_list]
threshold = 74.3  # Select pressure threshold

""" Todini index to evaluate water resilience

More info about Todini index:
<https://wntr.readthedocs.io/en/latest/apidoc/wntr.metrics.hydraulic.html?highlight=todini#wntr.metrics.hydraulic.todini_index> 

"""
todini = wntr.metrics.todini_index(head, pressure, demand, pump_flowrate, wn, threshold)


# Actual demand
actual_demand = results.node['demand'].loc[:, wn.junction_name_list]

# Actual demand sum for each node
actual_demand_tot = actual_demand.sum(axis=0)

# Expected demand
exp_demand = expected_demand(wn)

# Expected demand sum for each node
exp_demand_tot = exp_demand.sum(axis=0)

# Total demand all network
total_demand = exp_demand.values.sum()

# Exp tot - actual tot demand
sub = exp_demand_tot - actual_demand_tot
sub[sub < 0] = 0  # Evaluation adjustment

# Failure magnitude index
failure_magnitude = sub.div(total_demand)
failure_magnitude[failure_magnitude < 0] = 0

print("")
print("Failure magnitude:")
print(failure_magnitude)

# Plot water network based on the severity of the failure magnitude
wntr.graphics.plot_network(wn, node_attribute=failure_magnitude)  # Darker colors correspond to higher failure values
plt.title("Failure magnitude")


# CNA metrics plot and their correlation with failure magnitude

# Betweenness centrality
bce = pd.Series(betweenness_centrality)
wntr.graphics.plot_network(wn, node_attribute=betweenness_centrality)
plt.title("BCE")
print("")
print("Correlation Failure Magnitude/BCE:", round(failure_magnitude.corr(bce, method='pearson'), 3))

# Node degree
degree = pd.Series(node_degree)
wntr.graphics.plot_network(wn, node_attribute=degree)
plt.title("Node degree")
print("")
print("Correlation Failure Magnitude/node degree:", round(failure_magnitude.corr(degree, method='pearson'), 3))

# Closeness centrality
close_cen = pd.Series(closeness_centrality)
wntr.graphics.plot_network(wn, node_attribute=close_cen)
plt.title("Closeness centrality")
print("")
print("Correlation Failure Magnitude/Closeness Centrality:",
      round(failure_magnitude.corr(close_cen, method='pearson'), 3))

# Clustering coefficient
cluster_coeff = pd.Series(clustering_coeff)
wntr.graphics.plot_network(wn, node_attribute=cluster_coeff)
plt.title("Clustering coefficient")
print("")
print("Correlation Failure Magnitude/Clustering Coefficient:",
      round(failure_magnitude.corr(cluster_coeff, method='pearson'), 3))

# Eigenvector centrality
eig = pd.Series(eigenvector_centrality)
wntr.graphics.plot_network(wn, node_attribute=cluster_coeff)
plt.title("EIG centrality")
print("")
print("Correlation Failure magnitude/Eigenvector centrality:", round(failure_magnitude.corr(eig, method='pearson'), 3))

# Head water metrics and its correlation to CNA metrics and failure magnitude
head = np.mean(results.node['head'])
print("-----------------------------")
# print("Correlation head/Clustering Coefficient:",round(head.corr(cluster_coeff, method='pearson'),3))
# print("")
# print("Correlation head/Close_cen:",round(head.corr(close_cen, method='pearson'),3))
# print("")
# print("Correlation head/node_degree:",round(head.corr(degree, method='pearson'),3))
# print("")
# print("Correlation head/BCE:",round(head.corr(bce, method='pearson'),3))
# print("")
print("Correlation head/failure magnitude:", round(head.corr(failure_magnitude, method='pearson'), 3))
print("")
# print("Correlation head/Eigenvector centrality:",round(head.corr(eig, method='pearson'),3))
