#include <cmath>

#define MAX_CHILD 3
#define TREE_DEPTH 2
#define MAX_COPY 5
#define MAX_SUBTYPE ((((int)(pow(MAX_CHILD,TREE_DEPTH)))-1)/(MAX_CHILD-1))
#define NONLEAF ((((int)(pow(MAX_CHILD,TREE_DEPTH-1)))-1)/(MAX_CHILD-1))
#define MAX_H ((int)(pow(2,MAX_CHILD)))
