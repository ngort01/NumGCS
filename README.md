NumGCS
=======

NumGCS is a library that uses numerical optimization methods for Geometric Constraint Solving (GCS).
It is build on top of the [nlopt](https://github.com/stevengj/nlopt) library for numerical, nonlinear local and global optimization.


## Spatial Ontology
NumGCS implements a set of primitive spatial relations for points, lines and circles:
![](img/primitives.PNG?raw=true)

As described in [1], the primitive relations can be used to encode a large set of high-order spatial relations:
![](img/high_order.PNG?raw=true)





[1] Carl Schultz and Mehul Bhatt. A numerical optimisation based characterisation of spatial reasoning. In RuleML, 2016.
