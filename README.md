# RA_python (coming soon)
Python code for the Rectangular Aggregation (RA) technique used in [1707.05783](arxiv.org/abs/1707.05783) to mine the LHC dataset for hints of new physics. Based on the idea that, in the presence of a finely-binned parameter space, a signal of New Physics would populate multiple nearby bins.

It takes a *N*-dimensional parameter space over which exclusive bins are defined and generates aggregations of nearby bins into *super-bins* (or RA), and can then compute the (local) significance of each aggregation, which is a delta log-likelihood between the Standard Model hypothesis and the New Physics hypothesis where the signal only populates that one RA. 

Some examples based on 2017 CMS searches are provided.
