# matpower-simplification
A simpler, faster, and more readable version of matpower runpf; 
Deleting a lot of unnecessary steps to speed it up. CANNOT handle non-consecutive bus indices. 
CANNOT check Q limit and convert them into PQ buses.
In small cases it can be 10* faster than original; 
But when case gets bigger the difference will be smaller and smaller.
Requirements: matpower 7.0 installed, and the test and running is on MATLAB 2020b.
