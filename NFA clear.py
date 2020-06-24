from copy import deepcopy
from math import pow, sqrt
from operator import itemgetter
from sys import float_info
import matplotlib.pyplot as plt
from numpy import arccos, dot, pi, cross, asarray
from numpy.linalg import norm

### import points
def read_file_fence(dst, filePath):
    with open(filePath) as f:
        for line in f:
            dst.append([float(it) for it in line.split()])

### return true if line segments S1 and S2 intersect
# S1,S2 = segments [[x1,y1],[x2,y2]],[[x3,y3],[x4,y4]]
def segments_intersect(S1, S2):
    def orientation(p1,p2,p3):
        val = (p2[1]-p1[1])*(p3[0]-p2[0])-(p2[0]-p1[0])*(p3[1]-p2[1])
        return 0 if val == 0 else 1 if val > 0 else 2

    def on_segment(p1, p2, p3):
        return ((max(p1[0], p3[0]) >= p2[0] >= min(p1[0], p3[0])) and
                (max(p1[1], p3[1]) >= p2[1] >= min(p1[1], p3[1])))

    o1 = orientation(S1[0], S1[1], S2[0])
    o2 = orientation(S1[0], S1[1], S2[1])
    o3 = orientation(S2[0], S2[1], S1[0])
    o4 = orientation(S2[0], S2[1], S1[1])
    if (o1 != o2 and o3 != o4) or \
            (o1 == 0 and on_segment(S1[0], S2[0], S1[1])) or \
            (o2 == 0 and on_segment(S1[0], S2[1], S1[1])) or \
            (o3 == 0 and on_segment(S2[0], S1[0], S2[1])) or \
            (o4 == 0 and on_segment(S2[0], S1[1], S2[1])):
        return True
    return False

### return distance to segment A-B to point P
# S = segment index in figure, P = point index in points
def distance_to_segment(S, P):
    A = asarray(Points[Figure[S][0]])
    B = asarray(Points[Figure[S][1]])
    P = asarray(Points[P])
    """ segment line AB, point P, where each one is an array([x, y]) """
    if all(A == P) or all(B == P):
        return 0
    if arccos(dot((P - A) / norm(P - A), (B - A) / norm(B - A))) > pi / 2:
        return norm(P - A)
    if arccos(dot((P - B) / norm(P - B), (A - B) / norm(A - B))) > pi / 2:
        return norm(P - B)
    return norm(cross(A - B, A - P)) / norm(B - A)

### return distance between given 2 points as [x,y]
# p1, p2 = point index in points
def dist(p1, p2):
    p1 = Points[p1]
    p2 = Points[p2]
    return sqrt(pow(p1[0] - p2[0], 2) + pow(p1[1] - p2[1], 2))

### recalculate distances on given row (edge)
# put -1 on used points OR on points which cause intersection
def matrix_recalculate_row(row_1):
    for point in range(len(Points)):
        if Points_Visited[point]:
            distance_matrix[row_1][point] = -1
        else:
            check_intersection_with = (row_1 - 1) if row_1 != 0 else len(distance_matrix)-1
            if segments_intersect(
                    [Points[Figure[check_intersection_with][0]],Points[Figure[check_intersection_with][1]]],
                    [Points[Figure[row_1][1]],Points[point]]):
                distance_matrix[row_1][point] = -1
                continue
            check_intersection_with = (row_1 + 1) if row_1 != len(distance_matrix)-1 else 0
            if segments_intersect(
                    [Points[Figure[check_intersection_with][0]],Points[Figure[check_intersection_with][1]]],
                    [Points[Figure[row_1][0]],Points[point]]):
                distance_matrix[row_1][point] = -1
                continue

            distance_matrix[row_1][point] = distance_to_segment(row_1, point)

### cross out from matrix given point, by making distance -1 for each row for this point
def matrix_cross_out(point):
    for row_2 in range(len(distance_matrix)):
        distance_matrix[row_2][point] = -1

### return minimum distance in matrix distances as [point,edge]
# if there are few edges with equal distances, then answer is [[point_1,edge_1],[point_2,edge_2],...]
def matrix_get_minimum_distance():
    def sort_by_perimeter(case_1):
        segment_perimeters = [[dist(answer[case_1][0], Figure[answer[case_1][segment]][0]) +
                               dist(answer[case_1][0], Figure[answer[case_1][segment]][1]), answer[case][segment]]
                              for segment in range(1, len(answer[case_1]))]
        return [xxx[1] for xxx in sorted(segment_perimeters, key=itemgetter(0))]

    answer = []
    min = float_info.max

    for i_row in range(len(distance_matrix)):
        for i_point in range(len(Points)):
            if distance_matrix[i_row][i_point] == -1:
                continue
            else:
                local_tmp = distance_matrix[i_row][i_point]
                if local_tmp == min:
                    new_case = True
                    for case in range(len(answer)):
                        if i_point == answer[case][0]:
                            answer[case].append(i_row)
                            new_case = False
                            break
                    if new_case:
                        answer.append([i_point, i_row])
                elif local_tmp < min:
                    min = local_tmp
                    answer = [[i_point, i_row]]

    best_perimeter = float_info.max
    best_case = 0
    for case in range(len(answer)):
        potential = [sort_by_perimeter(case)[0]]
        potential.insert(0, answer[case][0])
        local_tmp = dist(potential[0], Figure[potential[1]][0]) + dist(potential[0], Figure[potential[1]][1])
        if local_tmp < best_perimeter:
            best_perimeter = local_tmp
            best_case = case
            answer[case] = potential
    return answer[best_case]

### Initialize first 2 edges of figure, called only once at the algorithm beginning.
# Current algorithm state: make init triangle from longest line and 1 random point
def set_init_figure(init_points):
    global Figure, Points, Points_Visited

    if len(init_points) < 2:
        print("Too few init points! Minimum 3.")
        exit(1)

    for it in range(len(init_points)):
        if Points_Visited[it]:
            print("Two exact init points found.")
            exit(1)
        else:
            Points_Visited[it] = True

        Figure.append([init_points[it], init_points[it + 1 if it + 1 != len(init_points) else 0]])
        distance_matrix.append([0] * len(Points))

    for it in range(len(init_points)):
        matrix_recalculate_row(it)

### Merge next point with figure
# args: LIST, args[0] is index of point, args[1] is index of figure segment
def pick_next(args):
    best_point = args[0]
    best_edge = args[1]

    Points_Visited[best_point] = True

    Figure.insert(best_edge+1,[best_point, Figure[best_edge][1]])
    Figure[best_edge][1] = best_point

    distance_matrix.insert(best_edge+1, [0]*len(Points))
    matrix_cross_out(best_point)
    matrix_recalculate_row(best_edge)
    matrix_recalculate_row(best_edge+1)

    print("{0:.1f}%".format(100 / len(Points) * len(Figure)))

### Calculate figure perimeter (TSP road length)
def get_tsp_length():
    length = 0.0
    for seg in range(len(Figure)):
        length += dist(Figure[seg][0],Figure[seg][1])
    return round(length, 1)

### 'Animate' point merges frame by frame
def animate(time = 0.1):
    mx = [Points[a[0]][0] for a in Figure]
    mx.append(mx[0])
    my = [Points[a[0]][1] for a in Figure]
    my.append(my[0])
    plt.cla()
    plt.plot(xx, yy, 'go')
    plt.plot(mx, my, 'ro', linestyle="solid", linewidth=2)
    plt.pause(time)
#######################################PROGRAM##########################################################################

Points, Figure, distance_matrix = [], [], []
read_file_fence(Points, 'berlin52.txt')
Points_Visited = [False] * len(Points)

xx = [a[0] for a in Points]
xx.append(xx[0])
yy = [a[1] for a in Points]
yy.append(yy[0])

set_init_figure([0,1,2])
while not all(Points_Visited):
    # I use there special dynamic matrix to store distances between segments and points.
    # You can use your own way to select and merge point to figure.
    pick_next(matrix_get_minimum_distance())
    animate(0.1)

print("Route: ", [p[0] for p in Figure])
print(len(Figure))
print(len(Points))
print("Length:", get_tsp_length())
plt.show()
