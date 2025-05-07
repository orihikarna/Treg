import builtins
import math

from . import kad, mat2

zero = (0, 0)


def sign(v):
    if v > 0:
        return +1
    if v < 0:
        return -1
    return 0

# vectors utility


def equal(a, b): return a[0] == b[0] and a[1] == b[1]
def length2(a):      ax, ay = a; return ax * ax + ay * ay
def length(a): return math.sqrt(length2(a))
def add(a, b):       ax, ay = a;   bx, by = b; return (ax + bx, ay + by)
def sub(a, b):       ax, ay = a;   bx, by = b; return (ax - bx, ay - by)
def dot(a, b):       ax, ay = a;   bx, by = b; return ax * bx + ay * by
def area(a, b):      ax, ay = a;   bx, by = b; return ax * by - ay * bx
def angle(a, b): return math.degrees(math.atan2(area(a, b), dot(a, b)))
def scale(scale, a, b=(0, 0)):   ax, ay = a;   bx, by = b; return (scale * ax + bx, scale * ay + by)


def distance2(a, b): return length2(sub(a, b))
def distance(a, b): return math.sqrt(distance2(a, b))
def proj(a, n): return scale(dot(a, n) / length2(n), n)
def perp(a, n): return sub(a, proj(a, n))
def round(a, ndigits=1): return (builtins.round(a[0], ndigits), builtins.round(a[1], ndigits))


def normalize(a):
    l = length(a)
    return scale(1 / l, a), l


def mult(m, a, b=(0, 0)):
    (ma, mc), (mb, md) = m
    ax, ay = a
    bx, by = b
    return (ma * ax + mb * ay + bx, mc * ax + md * ay + by)


def rotate(angle):
    th = math.radians(angle)
    co = math.cos(th)
    si = math.sin(th)
    return (co, si)


##
# Util
##

# lines
# k0 * t0 + p0 = k1 * t1 + p1
# t0 * k0 - t1 * k1 = p1 - p0
def find_intersection(p0, t0, p1, t1, tolerance=1):
    if True:  # exclude if nearly parallel
        degree = abs(angle(t0, t1))
        if degree < tolerance or 180 - tolerance < degree:
            return ((None, None), None, None)
        #print( 'degree = {}'.format( degree ) )
    A = (t0, t1)
    R = mat2.invert(A)
    B = sub(p1, p0)
    k0, k1 = mult(R, B)
    q0 = scale(k0, t0, p0)
    q1 = scale(-k1, t1, p1)
    dq = length(sub(q0, q1))
    if dq > 1e-3:
        print('ERROR |q0 - q1| = {}'.format(dq))
    return (q1, k0, -k1)


def combine_points(apnts, xpnt, bpnts):
    pnts = []
    for pnt in apnts:
        pnts.append(pnt)
    if xpnt is not None:
        pnts.append(xpnt)
    for pnt_b in reversed(bpnts):
        pnts.append(pnt_b)
    return pnts

##
# Interpolation
##


def interpolate_points_by_bezier(pnts, num_divs, debug):
    num_pnts = len(pnts)
    tmp = [None for _ in range(num_pnts)]
    curv = [pnts[0]]
    for i in range(1, num_divs):
        t = float(i) / num_divs
        s = 1 - t
        for n in range(num_pnts):
            tmp[n] = pnts[n]
        for L in range(num_pnts - 1, 0, -1):
            for n in range(L):
                tmp[n] = scale(s, tmp[n], scale(t, tmp[n+1]))
        curv.append(tmp[0])
    curv.append(pnts[-1])
    if debug:
        for pnt in pnts:
            kad.add_arc(pnt, add(pnt, (20, 0)), 360, 'F.Fab', 4)
    return curv


def h0(t): return (1 - t) * (1 - t) * (2 * t + 1)
def h1(t): return t * t * (3 - 2 * t)
def g0(t): return t * (1 - t) * (1 - t)
def g1(t): return t * t * (1 - t)


def interpolate_points_by_hermit_spline(apos, avec, bpos, bvec, num_divs, vec_scale):
    curv = [apos]
    for i in range(1, num_divs):
        t = float(i) / num_divs
        pos = add(scale(h0(t), apos), scale(h1(t), bpos))
        vec = add(scale(g0(t), (avec[1], -avec[0])), scale(g1(t), (bvec[1], -bvec[0])))
        pos = scale(vec_scale, vec, pos)
        curv.append(pos)
    curv.append(bpos)
    return curv

##
# Corner
##


def make_bezier_corner(corner, auvec, buvec, raidus, num_divs, debug):
    cef = (math.sqrt(0.5) - 0.5) / 3 * 8
    anchors = [
        scale(raidus,       auvec, corner),
        scale(raidus * cef, auvec, corner),
        scale(raidus * cef, buvec, corner),
        scale(raidus,       buvec, corner),
    ]
    return interpolate_points_by_bezier(anchors, num_divs, debug)


def make_arc_corner(corner, auvec, buvec, max_side_len, radius, num_divs, debug):
    if debug:
        kad.add_arc(corner, add(corner, (10, 0)), 360, 'F.Fab', 0.2)
        assert False
    theta = angle(auvec, buvec)
    # theta = builtins.round( theta * 10 ) / 10
    tan = math.tan((theta / 2) / 180 * math.pi)
    arc_radius = tan * max_side_len
    arc_radius = min(abs(arc_radius), radius) * sign(arc_radius)
    side_len = arc_radius / tan

    afoot = scale(side_len, auvec, corner)
    bfoot = scale(side_len, buvec, corner)
    aperp = (-auvec[1], auvec[0])

    if theta >= 0:
        arc_theta = 180 - theta
    else:
        arc_theta = -(180 + theta)
    arc_center = scale(arc_radius, aperp, afoot)
    # kad.add_arc( arc_center, add( arc_center, (0.6, 0) ), 360, 'F.Fab', 0.1 )

    pnts = [afoot]
    for i in range(1, num_divs):
        rot_mat = mat2.rotate(arc_theta / num_divs * i)
        pnt = scale(-arc_radius, mult(rot_mat, aperp), arc_center)
        pnts.append(pnt)
    pnts.append(bfoot)

    return pnts
