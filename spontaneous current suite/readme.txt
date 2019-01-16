%%%% READ ME FOR SPON. CURRENT SUITE %%%%%

-------- individual neuron analysis

1) Choose EPSC or IPSC. Assumes inward EPSC and outward IPSC.

2) Fill in parameters for recordings

3) Noise SD sets the threshold for detection above noise at "Noise SD" X the median absolute deviation of the noise using a 50ms rolling baseline (prior to event)

4) Template factor allows you to set a different "Noise SD" for putting together the template. Make this higher if you want less noise contaminating the template at the cost of subsampling from higher amplitude events to create the template.

5) Pick Folder and run

----- group comparisons

1) set parameters

2) experimental groups should be the names of the folders the data are stored in

3) Pick folder and run

-------

Not everything will run yet, still a work in progress. Individual neuron analysis should be completely functional and the group comparisons less so.
