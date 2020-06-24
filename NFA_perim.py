from math import pow, sqrt
from operator import itemgetter
from sys import float_info
from typing import List, Any

import matplotlib.pyplot as plt
from numpy import arccos, dot, pi, cross, asarray
from numpy.linalg import norm
from copy import deepcopy
from math import ceil, floor

### import points
def read_file_fence(filePath):
    global Fence_Points
    with open(filePath) as f:
        for line in f:
            Fence_Points.append([float(it) for it in line.split()])

def recalculate_qp():
    global quadrant_points
    quadrant_points = []
    for it in range(4):
        best = x.index(max(x)) if it == 0 else y.index(min(y)) if it == 1 else x.index(min(x)) if it == 2 else y.index(max(y))
        tmp_qp = []
        for pp in range(len(x)):
            if Fence_Points[pp][it % 2] == Fence_Points[best][it % 2]:
                tmp_qp.append(pp)

        if len(tmp_qp) < 2:
            quadrant_points.append(tmp_qp)
            continue

        tmp_d = {}
        for j in tmp_qp:
            tmp_d[j] = Fence_Points[j][it % 2 - 1]
        if it == 0:
            tmp_qp = [item[0] for item in sorted(tmp_d.items(), key = lambda kv:(kv[1], kv[0]))[::-1]]
        elif it == 1:
            tmp_qp = [item[0] for item in sorted(tmp_d.items(), key = lambda kv:(kv[1], kv[0]))][::-1]
        elif it == 2:
            tmp_qp = [item[0] for item in sorted(tmp_d.items(), key = lambda kv:(kv[1], kv[0]))]
        else:
            tmp_qp = [item[0] for item in sorted(tmp_d.items(), key = lambda kv:(kv[1], kv[0]))]
        quadrant_points.append(tmp_qp)

### find points which form minimum perimeter to fit all unvisited points
def get_next_perimeter():
    global x, y, Fence_Points

    if len(Fence_Points) <= 3:
        ret = [it for it in range(len(Fence_Points))]
        Fence_Points = []
        return ret

    x = [Fence_Points[it][0] for it in range(len(Fence_Points))]
    y = [Fence_Points[it][1] for it in range(len(Fence_Points))]

    recalculate_qp()
    current_quadrant = 0
    perimeter = []
    perimeter += quadrant_points[0]
    key_points = [it[0] for it in quadrant_points]

    while perimeter[0] != perimeter[-1] or len(perimeter) == 1:
        current_point = perimeter[-1]

        if current_quadrant == 0:
            possible_points = [it for it in range(len(x)) if y[it] < y[current_point] and x[it] < x[current_point]]
        elif current_quadrant == 1:
            possible_points = [it for it in range(len(x)) if y[it] > y[current_point] and x[it] < x[current_point]]
        elif current_quadrant == 2:
            possible_points = [it for it in range(len(x)) if y[it] > y[current_point] and x[it] > x[current_point]]
        else:
            possible_points = [it for it in range(len(x)) if y[it] < y[current_point] and x[it] > x[current_point]]

        if len(possible_points) == 0:
            if current_point == quadrant_points[current_quadrant+1][0]:
                perimeter += quadrant_points[current_quadrant+1][1:]
            current_quadrant += 1
            continue

        next_point = -1
        best_tg = float_info.max
        if current_quadrant % 2 == 0:
            for candidate in possible_points:
                if x[current_point]-x[candidate] == 0:
                    next_point = candidate
                    break
                if abs(x[current_point]-x[candidate])/abs(y[current_point]-y[candidate]) < best_tg:
                    next_point = candidate
                    best_tg = abs(x[current_point]-x[candidate])/abs(y[current_point]-y[candidate])
        else:
            for candidate in possible_points:
                if y[current_point]-y[candidate] == 0:
                    next_point = candidate
                    break
                if abs(y[current_point]-y[candidate])/abs(x[current_point]-x[candidate]) < best_tg:
                    next_point = candidate
                    best_tg = abs(y[current_point]-y[candidate])/abs(x[current_point]-x[candidate])

        if current_quadrant < 3 and next_point == key_points[current_quadrant+1]:
            perimeter += quadrant_points[current_quadrant+1]
        else:
            perimeter.append(next_point)
    for p in sorted(perimeter[:-1])[::-1]:
        Fence_Points.pop(p)
    return perimeter

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
    global Figure, Points

    for it in range(len(init_points)):
        Figure.append([init_points[it], init_points[it + 1 if it + 1 != len(init_points) else 0]])
        distance_matrix.append([0] * len(Points))

    for it in range(len(init_points)):
        matrix_recalculate_row(it)

### Pick up next point to add to figure with best new perimeter in case there are n>1 possible segments candidates
def pick_next():
    best_case = matrix_get_minimum_distance()
    best_point = best_case[0]
    best_edge = best_case[1]

    Points_Visited[best_point] = True

    Figure.insert(best_edge+1,[best_point, Figure[best_edge][1]])
    Figure[best_edge][1] = best_point

    distance_matrix.insert(best_edge+1, [0]*len(Points))
    matrix_cross_out(best_point)
    matrix_recalculate_row(best_edge)
    matrix_recalculate_row(best_edge+1)

    # print(Figure)
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
    plt.plot(perimp[0], perimp[1], 'bo', linestyle="None")
    plt.pause(time)

### Normalize points sequence
def translate_perimeter(my_perim):
    translated = []
    for per in my_perim:
        for pindex in range(len(Points)):
            if Points[pindex][0] == x[per] and Points[pindex][1] == y[per]:
                translated.append(pindex)
    return translated

def get_it(lis):
    if len(lis) <= 2:
        return lis
    return [lis[0],lis[-1],lis[floor(len(lis)/2)]]+get_it(lis[1:floor(len(lis)/2)])+get_it(lis[floor(len(lis)/2)+1:-1])
#######################################PROGRAM##########################################################################

quadrant_points, x, y, Fence_Points, Figure, distance_matrix, PL = [], [], [], [], [], [], []

read_file_fence('berlin52.txt')
Points = deepcopy(Fence_Points)

Points_Visited = [True] * len(Points)
xx = [a[0] for a in Points]
xx.append(xx[0])
yy = [a[1] for a in Points]
yy.append(yy[0])

##################################GENERATE_SEQUENCE_ZONE################################################################

perim = translate_perimeter(get_next_perimeter()[:-1])
while len(perim):
    PL.append(perim)
    perim = list(set(translate_perimeter(get_next_perimeter()[:-1])))

tmp, tmp1 = [], []
for i in range(ceil(len(PL)/2)):
    a = PL[i]
    b = PL[::-1][i]
    if len(PL) % 2 and i == ceil(len(PL)/2) - 1:
        tmp.append(a)
    else:
        tmp.append(a)
        tmp.append(b)
a = PL[::3]
b = PL[1::3]
c = PL[2::3]

for i in range(len(a)):
    tmp1.append(a[i])
    tmp1 += [c[i]] if i < len(c) else []
    tmp1 += [b[i]] if i < len(b) else []

# NFA combined with maximum perimeter idea with different strategies of choosing perimeters sequence
# Numbers show how perimeters were chosen

# BEST RESULTS I FOUND######## berlin52 # bier127  # tsp1000 #   tsp500  #   tsp250  #
# 123456789 -> 195 243 687
PL = get_it(PL)              #  7963.2  # 129885.4 # 27366.4  # 14469.8 #  14469.8  #
# 123456 -> 132465
# PL = tmp1                  #  8465.3  # 140433.6 # 30477.1  # 108792.2 #  14991.0  #
# 123456 -> 162534
# PL = tmp                   #  8361.8  # 138262.5 # 30260.3  # 107937.0 #  14994.4  #
# 123456 -> 135246
# PL = PL[::2] + PL[1:][::2] #  8522.6  # 138455.3 # 28642.2  # 101697.7 #  14454.9  #
# 123456 -> 654321
# PL = PL[::-1]              #  9550.0  # 149537.7 # 31842.3  # 112212.9 #  17209.4  #
# 123456 -> 123456
# PL = PL                    #  8749.1  # 142706.6 # 30751.8  # 108992.5 # 15016.7   #

#shuffle(PL)
################################END_GENERATE_SEQUENCE_ZONE##############################################################

set_init_figure(PL[0])

for perim in PL[1:]:
    perimp = [[Points[per][0] for per in perim],[Points[per][1] for per in perim]]
    for i in perim:
        Points_Visited[i] = False

    for i in range(len(distance_matrix)):
        matrix_recalculate_row(i)

    for iterator in range(len(perim)):
        pick_next()
    animate(0.1)

print("Route: ", [p[0] for p in Figure])
print(len(Figure))
print(len(Points))
print("Length:", get_tsp_length())
plt.show()
