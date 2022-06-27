# WDNsAnalysis
Failure analysis of Water Distribution Networks using WNTR and complex network theory methods

• Water distribution networks modelling based on Complex Network Theory

------------------------------------------------------------------------------------------------------

• About the project:

  The main goal of this project is the construction of a surrogate model of Water Distribution
  Network (WDN) applying Complex Network Theory principles. 
  This work combines concepts from water network analysis with elements from the analysis of complex networks and graphs.
  The proposed methodology is based on the Water Network Tool for Resilience (WNTR) tool, 
  which can generate a NetworkX data object that stores network connectivity as a graph. 
  The ability to easily integrate NetworkX with WNTR facilitates the use of numerous standard graph algorithms, 
  including algorithms that describe network structure.
  A NetworkX graph generated from a water network model stores the start and end node of each link, node coordinates, 
  and node and link types (i.e., tank, reservoir, valve). NetworkX includes numerous methods to analyze the structure of complex networks.

  The structural elements of water networks (e.g. valves, pumps, reservoirs, junctions and pipes) 
  are converted into the complex equivalent graph consisting of links and nodes, pipe length was used as graph weight.

  The project implements algorithms to evaluate the correlation between complex network theory metrics 
  and water distribution metrics.
  The calculations performed, evaluate the "small world" effect, graph theory measures for the 
  analysis of networks and hydraulic paramethers in relation to the intensity and duration of failures,
  for example, deficit in the fulfilled water demand, due to insufficient pressure in network. 

  The simulation paramethers, as well as graph metrics, can be tuned with values of structural elements like 
  nodes elevation, water flow, water quality, operational status, pressure, simulation time and others.

  Used failure indexes: 
  Failure duration: Node-related measure - Average duration of failure occurred to a specific node during simulation time

  Failure magnitude: Failure magnitude at nodes evaluates how bad the water shortage is, it is expressed
                     as the percentage of unmet water demand

  The results showed an inverse correlation between the "small world" property and the failure duration index. 

------------------------------------------------------------------------------------------------------

• Built with:

    Python (Numpy - Pandas - Matplotlib - scikit-learn - Scipy)
    Water Network Tool for Resilience (WNTR) (https://wntr.readthedocs.io/en/latest/index.html)
    NetworkX (https://networkx.org/)

------------------------------------------------------------------------------------------------------

• Project Folder:

  + WDNsAnalysis/networks: Folder containing benchmark networks
  + failure_duration.py: Script evaluating failure duration
  + falure_magnitude.py: Script evaluating failure magnitude
  + fail_correlation.py: Script evaluating correlation among metrics (main project file)
  + presentation.pdf: Project overview and results presentation

------------------------------------------------------------------------------------------------------

• Getting started:

   Prerequisites: Python, WNTR

------------------------------------------------------------------------------------------------------

• Installation:

1) Install WNTR tool in the local environment (https://wntr.readthedocs.io/en/latest/installation.html)

2) Choose a network from the "networks" folder or import an .inp file

3) Run the file with the selected network as input file





