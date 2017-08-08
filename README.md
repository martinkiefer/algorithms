# algorithms
Small Python-implementations of algorithms, implemented out of curiosity.

## Data Stream Processing
Algorithms do perform approximate or exact analyses over data streams

### Munro-Paterson's Algorithm (1978) for p-Pass Rank Selection 
Munro and Paterson were processing data streams before it was cool. Their publication is considered to be the first data stream processing publication and contains algorithms and upper bounds for sorting and rank selection over data streams / read-only tapes.

I implemented their "Multi-Pass Algorithm for Selection" (Section 3) which finds the element with a given rank in multiple passes over a data stream by using a deterministic sampling scheme (src/MunroPaterson.py). If the algorithm is not iterated until convergence, it can provide us with hard lower and upper bounds for the element with the given rank. The algorithm is invariant in terms of the magnitude of the used values as it is only based on pairwise comparisons and sorting. 

Reference:
J. I. Munro and M. S. Paterson. 1978. Selection and sorting with limited storage. In Proceedings of the 19th Annual Symposium on Foundations of Computer Science (SFCS '78). IEEE Computer Society, Washington, DC, USA, 253-258. DOI=http://dx.doi.org/10.1109/SFCS.1978.32
