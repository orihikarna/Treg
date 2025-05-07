import math
from . import kad, vec2

CIDX, POS, NIDX = range(3)


def connect_line_ends(line_ends, other_ends, curv_idx, pos, idx, layer, width):
    # overwrite new line_end
    line_ends[idx] = (curv_idx, pos, 0)

    # connect with up/down neighbors
    for nidx in [idx-1, idx+1]:
        if nidx < 0 or len(line_ends) <= nidx:
            continue
        # up / down neighbor
        neib = line_ends[nidx]
        if not neib:  # no neibor
            continue
        if (neib[NIDX] & (1 << idx)) != 0:  # already connected with this index
            #    L1-----R1      L1----
            #    C  ...---^^^ <-- cannot connect L1(right) and L2,
            #     L2--------R2    because L2 is already connected with L1(left)
            continue

        # The lines are drawn from left to right.
        # For the right line ends, reject the neighbor if there is a left end after the neighbor right end.
        # For the left line ends, reject the neighbor if there is a right end after the neighbor left end.
        if other_ends[nidx] and other_ends[nidx][CIDX] > neib[CIDX]:
            #    L1----------------------------R1
            #    C             .....-----^^^^^    <-- cannot connect R1 and R2(right),
            #     L2---------R2        L2---------R2  becuase L2(right) is on the right of R2(left)
            #          neib[CIDX] left[CIDX]
            continue
            pass
        # connect
        curr = line_ends[idx]
        line_ends[idx] = (curr[CIDX], curr[POS], curr[NIDX] | (1 << nidx))
        line_ends[nidx] = (neib[CIDX], neib[POS], neib[NIDX] | (1 << idx))
        kad.add_line(curr[POS], neib[POS], layer, width)


def load_distance_image(path):
    with open(path) as fin:
        data = fin.readlines()
        w = int(data[0])
        h = int(data[1])
        print(w, h)
        dist = []
        for y in range(h):
            vals = data[y+2].split(',')
            del vals[-1]
            if len(vals) != w:
                print('Error: len( vals )({}) != w({})'.format(len(vals), w))
                break
            vals = map(lambda v: float(v) / 100.0, vals)
            dist.append(vals)
    return (dist, w, h)


def get_distance(dist_image, pnt):
    dist, w, h = dist_image
    x, y = vec2.scale(10, vec2.sub(pnt, (Edge_CX, Edge_CY)))
    x, y = vec2.add((x, y), (w / 2.0, h / 2.0))
    x = min(max(x, 0.), w - 1.01)
    y = min(max(y, 0.), h - 1.01)
    nx = int()
    ny = int(math.floor(y))
    # if nx < 0 or w <= nx + 1 or ny < 0 or h <= ny + 1:
    #     return 0
    dx = x - nx
    dy = y - ny
    d = dist[ny][nx] * (1 - dx) * (1 - dy)
    d += dist[ny][nx+1] * dx * (1 - dy)
    d += dist[ny+1][nx] * (1 - dx) * dy
    d += dist[ny+1][nx+1] * dx * dy
    return d
