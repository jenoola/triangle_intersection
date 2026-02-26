This repository contains code and tests for a C++ implementation of detection of a triangle intersection in 3D space. The function *have_intersection* returns **true** if an intersection is detected and **false** otherwise.

Main algorithm is implemented according to *Olivier Devillers, Philippe Guigue. Faster Triangle-Triangle Intersection Tests. RR-4488, INRIA. 2002. inria-00072100*

Co-planar intersections are handled by impelementing *Möller–Trumbore ray-triangle intersection algorithm*. Degenerate cases (such as when a triangle turns out to be a point or a segment) are also handled.

