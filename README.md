# New Figure Algorithm V1.0
New greedy algorithm for symmetric (undirected) [TSP problem](https://en.wikipedia.org/wiki/Travelling_salesman_problem) is on development state.<br>
In the repository, you can find working code written in Python and some popular TSP maps.
#### You can use and modify my implementation if you want, but I'm going to explain the main idea, not the code. The code is a bit messed up and not written for explanation, because it's I continue improving it. Your implementation, of course, may be different.
## General Idea
The main idea is to extend the cycle by merging the nearest free vertexes to the existing cycle one by one.
Every iteration cycle is changing after a merge, so recalculation of distances between each segment of a cycle and free vertexes is needed.
The algorithm uses a special distance matrix that stores information about visited vertexes (vertexes which are already merged with graph) and
distances between each segment and free vertex. <br>
#### Example on berlin52 with three random init points:
![berlin52](https://imgur.com/H7sLjRC.gif)
### Figure array
To store a current path, the algorithm uses a 2D array ( [n][2] ) which stores every segment, forming a path. A segment is a 2 element subarray that stores indexes of two points, forming this segment.
### Distance Matrix
To reduce complexity, there is a special matrix that stores distances between figure segments and every point. Row index corresponds to the segment index in *Figure array* and the column index corresponds to the vertex index. <b> If a vertex is already present in figure or merging it with segment causes intersection, algorithm ignores this combination of vertex and segment by marking it with -1 value. Otherwise set value of the distance between segment and vertex. </b>
### Step by step
1. Initialize path with a minimum of 3 vertexes.
2. Recalculate distance matrix.
3. Check if any vertexes are free. If exist vertexes that are not in the path, then continue algorithm, else finish algorithm.
4. Find row and column with minimum non-negative value and merge selected vertex with the corresponding segment. Go to step 2.
#### Example graph
1.
![1](https://i.imgur.com/pY3lVba.jpg)
2.
![2](https://i.imgur.com/toPlWMQ.jpg)
3.
![3](https://i.imgur.com/3Oeq9M7.jpg)
4.
![4](https://i.imgur.com/w45qxQY.jpg)
#### Example on tsp250 with three random init points:
![TSP250](https://i.imgur.com/6GrcIRn.gif)

## Implementation of NFA with other ideas
This algorithm is flexible enough to combine it with other different ideas. <br>
For example, I have implemented an idea with a combination of perimeters with a minimum amount of vertexes which describe ( contain ) each other like tree rings. Then, I use NFA on vertexes from each circle in different combinations. The results are very different:
#### Merge sort - like sequence on TSP500
Length = 14469.8
![TSP500](https://i.imgur.com/lTqr4Qo.gif)

#### Regular from biggest to smallest circle sequence on TSP500
Length = 108992.5
![TSP500](https://i.imgur.com/COcxGYW.gif)

# Other examples:
#### berlin52 with circles implementation:
![berlin52](https://i.imgur.com/Maw9co6.gif)

# TODO:
- Consider adding a 3-opt
- Keep existing edges and distances of unjoined nodes to edges in ordered lists. Replace matrix with ordered list.
