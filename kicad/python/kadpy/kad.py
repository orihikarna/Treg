import math
from . import pnt, vec2, mat2
import pcbnew

UnitMM = True
PointDigits = 3
inf = 1e9


def scalar_to_unit(v, mm_or_mils):
    if mm_or_mils:
        return pcbnew.FromMM(v)
    else:
        return pcbnew.FromMils(v)


def to_ANGLE(angle):
    if type(angle) == pcbnew.EDA_ANGLE:
        return angle
    else:
        return pcbnew.EDA_ANGLE(angle)


pcb = pcbnew.GetBoard()

# wire types
Straight = 0
Directed = 1
ZigZag = 2


# draw corner types
Line = 0
Linear = 1
Bezier = 2
BezierRound = 3
Round = 4
LinearRound = 5
Spline = 6


##
# Drawings
##
def add_line_rawunit(a, b, layer="Edge.Cuts", width=2):
    if a[0] == b[0] and a[1] == b[1]:
        print("add_line_rawunit: identical", a)
        return None
    line = pcbnew.PCB_SHAPE(pcb)
    line.SetShape(pcbnew.SHAPE_T_SEGMENT)
    line.SetStart(a)
    line.SetEnd(b)
    line.SetLayer(pcb.GetLayerID(layer))
    line.SetWidth(pcbnew.FromMils(width))
    pcb.Add(line)
    return line


def add_line(a, b, layer="Edge.Cuts", width=2):
    pnt_a = pnt.to_unit(vec2.round(a, PointDigits), UnitMM)
    pnt_b = pnt.to_unit(vec2.round(b, PointDigits), UnitMM)
    if True:
        aa = pnt.from_unit(pnt_a, UnitMM)
        bb = pnt.from_unit(pnt_b, UnitMM)
        length = vec2.distance(bb, aa)
        if length == 0:
            print("add_line: identical", a)
            return
        elif length < 0.001:
            # print( length, a, b )
            None
    line = pcbnew.PCB_SHAPE(pcb)
    line.SetShape(pcbnew.SHAPE_T_SEGMENT)
    line.SetStart(pnt.to_VEC2I(pnt_a))
    line.SetEnd(pnt.to_VEC2I(pnt_b))
    line.SetLayer(pcb.GetLayerID(layer))
    line.SetWidth(scalar_to_unit(width, UnitMM))
    pcb.Add(line)
    return line


def add_lines(curv, layer, width):
    for idx, pnt in enumerate(curv):
        if idx == 0:
            continue
        add_line(curv[idx - 1], pnt, layer, width)


def add_arc(ctr, pos, angle, layer="Edge.Cuts", width=2):
    pnt_ctr = pnt.to_unit(vec2.round(ctr, PointDigits), UnitMM)
    pnt_pos = pnt.to_unit(vec2.round(pos, PointDigits), UnitMM)
    arc = pcbnew.PCB_SHAPE(pcb)
    arc.SetShape(pcbnew.SHAPE_T_ARC)
    arc.SetCenter(pnt_ctr)
    arc.SetStart(pnt_pos)
    arc.SetArcAngleAndEnd(10 * angle, True)
    arc.SetLayer(pcb.GetLayerID(layer))
    arc.SetWidth(scalar_to_unit(width, UnitMM))
    pcb.Add(arc)
    return arc


def add_arc2(ctr, pos, end, angle, layer="Edge.Cuts", width=2):
    pnt_ctr = pnt.to_unit(vec2.round(ctr, PointDigits), UnitMM)
    pnt_pos = pnt.to_unit(vec2.round(pos, PointDigits), UnitMM)
    pnt_end = pnt.to_unit(vec2.round(end, PointDigits), UnitMM)
    arc = pcbnew.PCB_SHAPE(pcb)
    arc.SetShape(pcbnew.SHAPE_T_ARC)
    arc.SetCenter(pnt.to_VEC2I(pnt_ctr))
    arc.SetStart(pnt.to_VEC2I(pnt_pos))
    arc.SetArcAngleAndEnd(to_ANGLE(angle), True)
    arc.SetLayer(pcb.GetLayerID(layer))
    arc.SetWidth(scalar_to_unit(width, UnitMM))
    pcb.Add(arc)
    arc_start = arc.GetStart()
    arc_end = arc.GetEnd()
    if arc_start[0] != pnt_pos[0] or arc_start[1] != pnt_pos[1]:
        print("add_arc2: arc_start != pnt_pos !!")
    if arc_end[0] != pnt_end[0] or arc_end[1] != pnt_end[1]:
        # print( 'arc_end != pnt_end !!' )
        add_line_rawunit(arc_end, pnt_end, layer, width)
    return arc


def add_arc3(ctr, pos, end, angle, layer="Edge.Cuts", width=2):
    vec = vec2.sub(pos, ctr)
    radius = vec2.length(vec)
    angle0 = vec2.angle((1, 0), vec)
    num_angle = int(math.ceil(angle / 3.0))  # degrees
    num_length = int(math.ceil(radius * math.radians(angle) / 1.0))  # mm
    num = max(num_angle, num_length)
    dangle = angle / num
    curv = [pos]
    for n in range(1, num):
        pnt = vec2.scale(radius, vec2.rotate(angle0 + dangle * n), ctr)
        curv.append(pnt)
    curv.append(end)
    add_lines(curv, layer, width)
    return curv


def add_text(
    pos,
    angle,
    string,
    layer="F.SilkS",
    size=(1, 1),
    thick=5,
    hjustify=None,
    vjustify=None,
):
    text = pcbnew.PCB_TEXT(pcb)
    text.SetPosition(pnt.to_VEC2I(pnt.to_unit(vec2.round(pos, PointDigits), UnitMM)))
    text.SetTextAngle(to_ANGLE((angle)))
    text.SetText(string)
    text.SetLayer(pcb.GetLayerID(layer))
    text.SetTextSize(pnt.to_VEC2I(pcbnew.wxSizeMM(size[0], size[1])))
    text.SetTextThickness(scalar_to_unit(thick, UnitMM))
    text.SetMirrored(layer[0] == "B")
    if hjustify != None:
        text.SetHorizJustify(hjustify)
    if vjustify != None:
        text.SetVertJustify(vjustify)
    pcb.Add(text)
    return text


##
## Tracks & Vias
##


def add_track(a, b, net, layer, width):
    pnt_a = pnt.to_VEC2I(pnt.to_unit(vec2.round(a, PointDigits), UnitMM))
    pnt_b = pnt.to_VEC2I(pnt.to_unit(vec2.round(b, PointDigits), UnitMM))
    track = pcbnew.PCB_TRACK(pcb)
    track.SetStart(pnt_a)
    track.SetEnd(pnt_b)
    if net != None:
        track.SetNet(net)
    track.SetLayer(layer)
    track.SetWidth(scalar_to_unit(width, UnitMM))
    track.SetLocked(True)
    pcb.Add(track)
    return track


def add_via(pos, net, size):  # size [mm]
    pnt_ = pnt.to_VEC2I(pnt.to_unit(vec2.round(pos, PointDigits), UnitMM))
    via = pcbnew.PCB_VIA(pcb)
    via.SetPosition(pnt_)
    via.SetWidth(pcbnew.FromMM(size[0]))
    via.SetDrill(pcbnew.FromMM(size[1]))
    via.SetNet(net)
    # print( net.GetNetname() )
    via.SetLayerPair(pcb.GetLayerID("F.Cu"), pcb.GetLayerID("B.Cu"))
    via.SetLocked(True)
    pcb.Add(via)
    return via


# def add_via_on_pad( mod_name, pad_name, via_size ):
#     pos, net = get_pad_pos_net( mod_name, pad_name )
#     add_via( pos, net, via_size )


def add_via_relative(mod_name, pad_name, offset_vec, size_via):
    pos_via = calc_pos_from_pad(mod_name, pad_name, offset_vec)
    net = get_pad_net(mod_name, pad_name)
    return add_via(pos_via, net, size_via)


def get_via_pos(via):
    return pnt.from_unit(via.GetPosition(), UnitMM)


def get_via_pos_net(via):
    return pnt.from_unit(via.GetPosition(), UnitMM), via.GetNet()


##
# Module
##
def _get_mod_pos(mod):
    return pnt.from_unit(mod.GetPosition(), UnitMM)


def _get_mod_angle(mod):
    return mod.GetOrientation().AsDegrees()


def _get_mod_layer(mod):
    return mod.GetLayerName()


#
def get_mod(mod_name):
    return pcb.FindFootprintByReference(mod_name)


def get_mod_pos(mod_name):
    mod = get_mod(mod_name)
    return _get_mod_pos(mod)


def get_mod_angle(mod_name):
    mod = get_mod(mod_name)
    return _get_mod_angle(mod)


def get_mod_layer(mod_name):
    mod = get_mod(mod_name)
    return _get_mod_layer(mod)


def get_mod_pos_angle(mod_name):
    mod = get_mod(mod_name)
    return _get_mod_pos(mod), _get_mod_angle(mod)


def set_mod_pos_angle(mod_name, pos, angle):
    mod = get_mod(mod_name)
    if pos is not None:
        mod.SetPosition(pnt.to_VEC2I(pnt.to_unit(vec2.round(pos, PointDigits), UnitMM)))
    mod.SetOrientation(to_ANGLE(angle))
    return mod


# mods


def move_mods(base_pos, base_angle, mods):
    mrot = mat2.rotate(base_angle)
    for vals in mods:
        name, pos, angle = vals[:3]
        npos = vec2.add(base_pos, vec2.mult(mrot, pos))
        nangle = base_angle + angle
        if name is None:
            move_mods(npos, nangle, vals[3])
        else:
            set_mod_pos_angle(name, npos, nangle)


# pad


def _get_pad(mod, pad_name):
    return mod.FindPadByNumber(pad_name)


def _get_pad_pos(pad):
    return pnt.from_unit(pad.GetPosition(), UnitMM)


def _get_pad_net(pad):
    return pad.GetNet()


def _get_pad_pos_net_angle_layer(mod, pad):
    return (
        _get_pad_pos(pad),
        _get_pad_net(pad),
        _get_mod_angle(mod),
        _get_mod_layer(mod),
    )


#
def get_mod_pad(mod_name, pad_name):
    mod = get_mod(mod_name)
    pad = _get_pad(mod, pad_name)
    return mod, pad


def get_pad(mod_name, pad_name):
    mod = get_mod(mod_name)
    pad = _get_pad(mod, pad_name)
    return pad


def get_pad_pos(mod_name, pad_name):
    pad = get_pad(mod_name, pad_name)
    return _get_pad_pos(pad)


def get_pad_net(mod_name, pad_name):
    pad = get_pad(mod_name, pad_name)
    return _get_pad_net(pad)


def get_pad_pos_net_angle_layer(mod_name, pad_name):
    mod, pad = get_mod_pad(mod_name, pad_name)
    return _get_pad_pos_net_angle_layer(mod, pad)


def calc_relative_vec(mod_name, vec, pos=(0, 0)):
    mod = get_mod(mod_name)
    angle = _get_mod_angle(mod)
    layer = _get_mod_layer(mod)
    if layer == "B.Cu":
        vec = (vec[0], -vec[1])
    vec_relative = vec2.mult(mat2.rotate(angle), vec, pos)
    return vec_relative


def calc_pos_from_pad(mod_name, pad_name, offset_vec):
    pos_pad = get_pad_pos(mod_name, pad_name)
    pos_relative = calc_relative_vec(mod_name, offset_vec, pos_pad)
    return pos_relative


def calc_pos_from_mod(mod_name, offset_vec):
    pos_mod = get_mod_pos(mod_name)
    pos_relative = calc_relative_vec(mod_name, offset_vec, pos_mod)
    return pos_relative


##
# Wires
##


def add_wire_straight(pnts, net, layer, width, radii):
    # assert radius >= 0
    num_pnts = len(pnts)
    rpnts = []
    for idx, curr in enumerate(pnts):
        if idx == 0 or idx == num_pnts - 1:  # first or last
            rpnts.append(curr)
            continue
        radius = radii[idx]
        prev = pnts[idx - 1]
        next = pnts[idx + 1]
        avec = vec2.sub(prev, curr)
        bvec = vec2.sub(next, curr)
        alen = vec2.length(avec)
        blen = vec2.length(bvec)
        side_len = min(
            alen / 2 if idx - 1 > 0 else alen,
            blen / 2 if idx + 1 < num_pnts - 1 else blen,
        )
        if side_len < 10 ** (-PointDigits):
            rpnts.append(curr)
        else:
            num_divs = 15
            debug = False
            auvec = vec2.scale(1 / alen, avec)
            buvec = vec2.scale(1 / blen, bvec)
            # rpnts += vec2.make_bezier_corner( curr, auvec, buvec, length, num_divs, debug )
            rpnts += vec2.make_arc_corner(curr, auvec, buvec, side_len, radius, num_divs, debug)
    for idx, curr in enumerate(rpnts):
        if idx == 0:
            prev = rpnts[0]
            continue
        length = vec2.distance(prev, curr)
        if length > 0.01:
            add_track(prev, curr, net, pcb.GetLayerID(layer), width)
            prev = curr


# params: pos, (offset length, offset angle / arc center) x n, direction angle


def add_wire_offsets_directed(prms_a, prms_b, net, layer, width, radius, arc_ctr_mid=None):
    def _make_points_from_offsets(start_pos, offsets, angle, radius):
        pos = start_pos
        pnts = [pos]
        rads = [None]
        for n, (off_angle, off_len_or_arcctr) in enumerate(offsets):
            off_len = None
            off_udir = vec2.rotate(-off_angle)
            if type(off_len_or_arcctr) == type(()):  # tuple
                arc_ctr = off_len_or_arcctr
                vec = vec2.sub(arc_ctr, pos)
                arc_rad = vec2.length(vec2.perp(vec, off_udir))
                off_len = vec2.length(vec2.proj(vec, off_udir))
                if vec2.dot(vec, off_udir) < 0:
                    off_len *= -1
                #
                next_off_angle = offsets[n + 1][0] if n + 1 < len(offsets) else angle
                ctr_angle = off_angle - next_off_angle
                while ctr_angle > +180:
                    ctr_angle -= 360
                while ctr_angle < -180:
                    ctr_angle += 360
                ctr_angle = abs(ctr_angle) / 2
                #
                off_len += arc_rad * math.tan(ctr_angle / 180 * math.pi)
                if True:  # debug
                    # if vec2.distance( arc_ctr, (183, 133) ) < 10:# debug
                    # arc_rad = 0
                    # print( f'{arc_rad = :.2f}, {off_len = :.2f}' )
                    # print( f'{ctr_angle = :.2f}, {off_angle = :.1f}, {next_off_angle = :.1f}' )
                    # print( f'{pos = }, {arc_ctr = }' )
                    add_arc(arc_ctr, vec2.add(arc_ctr, (0.4, 0)), 360, "F.Fab", 0.2)
            else:
                off_len = off_len_or_arcctr
                arc_rad = radius
            pos = vec2.scale(off_len, off_udir, pos)
            pnts.append(pos)
            rads.append(arc_rad)
        return pnts, rads

    #
    pos_a, offsets_a, angle_a = prms_a
    pos_b, offsets_b, angle_b = prms_b
    pnts_a, rads_a = _make_points_from_offsets(pos_a, offsets_a, angle_a, radius)
    pnts_b, rads_b = _make_points_from_offsets(pos_b, offsets_b, angle_b, radius)
    #
    end_a = pnts_a[-1]
    end_b = pnts_b[-1]
    dir_a = vec2.rotate(-angle_a)
    dir_b = vec2.rotate(-angle_b)
    xpos, _, _ = vec2.find_intersection(end_a, dir_a, end_b, dir_b)
    if xpos[0] is None:
        print(f"{xpos = }, {angle_a = }, {angle_b = }, {end_a = }, {end_b = }")
    if arc_ctr_mid:
        add_arc(arc_ctr_mid, vec2.add(arc_ctr_mid, (0.4, 0)), 360, "B.Fab", 0.2)
        dist_a = vec2.length(vec2.perp(vec2.sub(arc_ctr_mid, end_a), dir_a))
        dist_b = vec2.length(vec2.perp(vec2.sub(arc_ctr_mid, end_b), dir_b))
        radius_mid = min(dist_a, dist_b)
        # print( f'{radius_mid = }' )
    else:
        radius_mid = radius
    #
    pnts = vec2.combine_points(pnts_a, xpos, pnts_b)
    rads = vec2.combine_points(rads_a, radius_mid, rads_b)
    add_wire_straight(pnts, net, layer, width, rads)


# params: parallel lines direction angle


def add_wire_zigzag(pos_a, pos_b, angle, delta_angle, net, layer, width, radius):
    mid_pos = vec2.scale(0.5, vec2.add(pos_a, pos_b))
    dir = vec2.rotate(-angle)
    mid_dir1 = vec2.rotate(-angle + delta_angle)
    mid_dir2 = vec2.rotate(-angle - delta_angle)
    _, ka1, _ = vec2.find_intersection(pos_a, dir, mid_pos, mid_dir1)
    _, ka2, _ = vec2.find_intersection(pos_a, dir, mid_pos, mid_dir2)
    mid_angle = (angle - delta_angle) if abs(ka1) < abs(ka2) else (angle + delta_angle)
    add_wire_offsets_directed((pos_a, [], angle), (mid_pos, [], mid_angle), net, layer, width, radius)
    add_wire_offsets_directed((pos_b, [], angle), (mid_pos, [], mid_angle), net, layer, width, radius)


def __wire_mod_sub(pos_a, angle_a, sign_a, pos_b, angle_b, sign_b, net, layer, width, prms):
    def _proc_directed_params(prms, angle, sign):
        offsets = []
        if type(prms) == type([]):  # array
            for off_angle, off_len_or_arcctr in prms[0:-1]:
                offsets.append((angle + off_angle * sign, off_len_or_arcctr))
            dir_angle = prms[-1] * sign
        else:
            dir_angle = prms * sign
        return offsets, dir_angle

    #
    if type(prms) == type(Straight) and prms == Straight:
        add_wire_straight([pos_a, pos_b], net, layer, width, None)
    elif prms[0] == Directed:
        prms_a, prms_b = prms[1:3]
        offsets_a, dir_angle_a = _proc_directed_params(prms_a, angle_a, sign_a)
        offsets_b, dir_angle_b = _proc_directed_params(prms_b, angle_b, sign_b)
        prms2_a = (pos_a, offsets_a, angle_a + dir_angle_a)
        prms2_b = (pos_b, offsets_b, angle_b + dir_angle_b)
        #
        radius = prms[3] if len(prms) > 3 else inf
        arc_ctr_mid = prms[4] if len(prms) > 4 else None
        add_wire_offsets_directed(prms2_a, prms2_b, net, layer, width, radius, arc_ctr_mid)
    elif prms[0] == ZigZag:
        dangle, delta_angle = prms[1:3]
        radius = prms[3] if len(prms) > 3 else inf
        add_wire_zigzag(
            pos_a,
            pos_b,
            angle_a + dangle * sign_a,
            delta_angle,
            net,
            layer,
            width,
            radius,
        )


def wire_mod_pads(tracks):
    def _get_pad_props(pad, mod):
        if type(pad) is pcbnew.PCB_VIA:  # pad is Via
            pos, net = get_via_pos_net(pad)
            if mod is not None:  # ref = mod
                angle = get_mod_angle(mod)
                layer = get_mod_layer(mod)
            else:
                angle = None
                layer = None
        else:  # pab is Pad
            pos, net, angle, layer = get_pad_pos_net_angle_layer(mod, pad)
        return pos, angle, layer, net

    #
    layer_FCu = "F.Cu"
    layer_BCu = "B.Cu"
    for track in tracks:
        if track == None:
            continue
        mod_a, pad_a, mod_b, pad_b, width, prms = track[:6]
        # check a & b
        pos_a, angle_a, layer_a, net_a = _get_pad_props(pad_a, mod_a)
        pos_b, angle_b, layer_b, net_b = _get_pad_props(pad_b, mod_b)
        # process None
        if angle_a == None:
            angle_a = angle_b
        if angle_b == None:
            angle_b = angle_a
        sign_a = +1 if layer_a == None or layer_a == layer_FCu else -1
        sign_b = +1 if layer_b == None or layer_b == layer_FCu else -1
        net = net_b if net_a == None else net_a
        layer = layer_b if layer_a == None else layer_a
        # validation
        if net == None:
            print("net is None")
        if layer == None:
            print("layer is None")
        if len(track) > 6:
            if track[6] == "Opp":
                layer = layer_BCu if layer == layer_FCu else layer_FCu
            else:
                layer = track[6]
        __wire_mod_sub(pos_a, angle_a, sign_a, pos_b, angle_b, sign_b, net, layer, width, prms)


##
# Drawwings
##
def removeDrawings():
    # pcb.DeleteZONEOutlines()
    for zone in reversed(pcb.Zones()):
        # print(zone)
        pcb.Remove(zone)
    for draw in pcb.GetDrawings():
        pcb.Delete(draw)


def removeTracksAndVias():
    # Tracks & Vias
    for track in pcb.GetTracks():
        delete = False
        if track.IsLocked():  # locked == placed by python
            delete = True
        # elif type( track ) is pcbnew.PCB_VIA:
        #     delete = True
        # elif type( track ) is pcbnew.PCB_TRACK:
        #     delete = True
        if delete:
            pcb.Delete(track)


def drawRect(pnts, layer, R=75, width=2):
    if R > 0:
        arcs = []
        for idx, a in enumerate(pnts):
            b = pnts[idx - 1]
            vec = vec2.sub(b, a)
            length = vec2.length(vec)
            delta = vec2.scale(R / length, vec)
            na = vec2.add(a, delta)
            nb = vec2.sub(b, delta)
            ctr = vec2.add(nb, (delta[1], -delta[0]))
            arcs.append((ctr, na, nb))
            add_line(na, nb, layer, width)
        for idx, (ctr, na, nb) in enumerate(arcs):
            arc = add_arc2(ctr, nb, arcs[idx - 1][1], -90, layer, width)
            # print( "na   ", pnt.mils2unit( vec2.round( na ) ) )
            # print( "start", arc.GetArcStart() )
            # print( "nb   ", pnt.mils2unit( vec2.round( nb ) ) )
            # print( "end  ", arc.GetArcEnd() )
    else:
        for idx, a in enumerate(pnts):
            b = pnts[idx - 1]
            add_line(a, b, layer, width)


# edge.cut drawing
def is_supported_round_angle(angle):
    # integer?
    if angle - round(angle) != 0:
        return False
    # multiple of 90?
    deg = int(angle)
    if (deg % 90) != 0:
        return False
    return True


def calc_bezier_corner_points(apos, avec, bpos, bvec, pitch=1, ratio=0.7):
    _, alen, blen = vec2.find_intersection(apos, avec, bpos, bvec)
    debug = False
    if alen <= 0:
        print("BezierCorner: alen = {} < 0, at {}".format(alen, apos))
        debug = True
    if blen <= 0:
        print("BezierCorner: blen = {} < 0, at {}".format(blen, bpos))
        debug = True
    # if debug:
    #     return []
    actrl = vec2.scale(alen * ratio, avec, apos)
    bctrl = vec2.scale(blen * ratio, bvec, bpos)
    ndivs = int(round((alen + blen) / pitch))
    pnts = [apos, actrl, bctrl, bpos]
    curv = vec2.interpolate_points_by_bezier(pnts, ndivs, debug)
    return curv


def calc_bezier_round_points(apos, avec, bpos, bvec, radius):
    _, alen, blen = vec2.find_intersection(apos, avec, bpos, bvec)
    debug = False
    if alen <= radius:
        print("BezierRound: alen < radius, {} < {}, at {}".format(alen, radius, apos))
        debug = True
    if blen <= radius:
        print("BezierRound: blen < radius, {} < {}, at {}".format(blen, radius, bpos))
        debug = True
    amid = vec2.scale(alen - radius, avec, apos)
    bmid = vec2.scale(blen - radius, bvec, bpos)
    # add_line( apos, amid, layer, width )
    # add_line( bpos, bmid, layer, width )
    angle = vec2.angle(avec, bvec)
    if angle < 0:
        # angle += 360
        angle *= -1
    ndivs = int(round(angle / 4.5))
    # print( 'BezierRound: angle = {}, ndivs = {}'.format( angle, ndivs ) )
    coeff = (math.sqrt(0.5) - 0.5) / 3 * 8
    # print( 'coeff = {}'.format( coeff ) )
    actrl = vec2.scale(alen + (coeff - 1) * radius, avec, apos)
    bctrl = vec2.scale(blen + (coeff - 1) * radius, bvec, bpos)
    pnts = [amid, actrl, bctrl, bmid]
    curv = [apos]
    for pos in vec2.interpolate_points_by_bezier(pnts, ndivs, debug):
        curv.append(pos)
    curv.append(bpos)
    return curv


def draw_corner(cnr_type, a, cnr_data, b, layer, width):
    apos, aangle = a
    bpos, bangle = b
    avec = vec2.rotate(aangle)
    bvec = vec2.rotate(bangle + 180)
    curv = None
    if cnr_type != Bezier and abs(vec2.dot(avec, bvec)) > 0.999:
        cnr_type = Line
        # print( avec, bvec )
    if cnr_type == Line:
        add_line(apos, bpos, layer, width)
    elif cnr_type == Linear:
        xpos, alen, blen = vec2.find_intersection(apos, avec, bpos, bvec)
        if cnr_data == None:
            add_line(apos, xpos, layer, width)
            add_line(bpos, xpos, layer, width)
        else:
            delta = cnr_data[0]
            # print( delta, alen, xpos )
            amid = vec2.scale(alen - delta, avec, apos)
            bmid = vec2.scale(blen - delta, bvec, bpos)
            add_line(apos, amid, layer, width)
            add_line(bpos, bmid, layer, width)
            add_line(amid, bmid, layer, width)
    elif cnr_type == Bezier:
        num_data = len(cnr_data)
        alen = cnr_data[0]
        blen = cnr_data[-2]
        ndivs = cnr_data[-1]
        apos2 = vec2.scale(alen, avec, apos)
        bpos2 = vec2.scale(blen, bvec, bpos)
        pnts = [apos, apos2]
        if num_data > 3:
            for pt in cnr_data[1 : num_data - 2]:
                pnts.append(pt)
        pnts.append(bpos2)
        pnts.append(bpos)
        curv = vec2.interpolate_points_by_bezier(pnts, ndivs)
        add_lines(curv, layer, width)
    elif cnr_type == BezierRound:
        radius = cnr_data[0]
        curv = calc_bezier_round_points(apos, avec, bpos, bvec, radius)
        add_lines(curv, layer, width)
    elif cnr_type in [Round, LinearRound]:
        radius = cnr_data[0]
        # print( 'Round: radius = {}'.format( radius ) )
        # print( apos, avec, bpos, bvec )
        xpos, alen, blen = vec2.find_intersection(apos, avec, bpos, bvec)
        # print( xpos, alen, blen )
        debug = False
        if not is_supported_round_angle(aangle):
            pass
            # print( 'Round: warning aangle = {}'.format( aangle ) )
            # debug = True
        if not is_supported_round_angle(bangle):
            pass
            # print( 'Round: warning bangle = {}'.format( bangle ) )
            # debug = True
        angle = vec2.angle(avec, bvec)
        angle = round(angle * 10) / 10
        tangent = math.tan(math.radians(abs(angle) / 2))
        side_len = radius / tangent
        if alen < side_len:
            print("Round: alen < side_len, {} < {}".format(alen, side_len))
            debug = True
        if blen < side_len:
            print("Round: blen < side_len, {} < {}".format(blen, side_len))
            debug = True
        if debug:
            add_arc(xpos, vec2.add(xpos, (10, 0)), 360, layer, width)
            return b, curv
        # print( 'angle = {}, radius = {}, side_len = {}, tangent = {}'.format( angle, radius, side_len, tangent ) )

        amid = vec2.scale(-side_len, avec, xpos)
        bmid = vec2.scale(-side_len, bvec, xpos)
        add_line(apos, amid, layer, width)
        add_line(bpos, bmid, layer, width)

        aperp = (-avec[1], avec[0])
        if angle >= 0:
            ctr = vec2.scale(-radius, aperp, amid)
            if cnr_type == Round:
                add_arc2(ctr, bmid, amid, 180 - angle, layer, width)
            elif cnr_type == LinearRound:
                add_arc3(ctr, bmid, amid, 180 - angle, layer, width)
        else:
            ctr = vec2.scale(+radius, aperp, amid)
            if cnr_type == Round:
                add_arc2(ctr, amid, bmid, 180 + angle, layer, width)
            elif cnr_type == LinearRound:
                add_arc3(ctr, amid, bmid, 180 + angle, layer, width)
    elif cnr_type == Spline:
        vec_scale = cnr_data[0]
        ndivs = int(round(vec2.distance(apos, bpos) / 2.0))
        auvec = vec2.rotate(aangle + 90)
        buvec = vec2.rotate(bangle - 90)
        curv = vec2.interpolate_points_by_hermit_spline(apos, auvec, bpos, buvec, ndivs, vec_scale)
        add_lines(curv, layer, width)
    return b, curv


def draw_closed_corners(corners, layer, width):
    curvs = []
    a = corners[-1][0]
    for b, cnr_type, cnr_data in corners:
        a, curv = draw_corner(cnr_type, a, cnr_data, b, layer, width)
        curvs.append(curv)
    return curvs


# zones


def _add_area(pnts, layer, net_name):
    net = -1 if net_name is None else pcb.FindNet(net_name).GetNetCode()
    pnts = [pnt.to_unit(vec2.round(pt, PointDigits), UnitMM) for pt in pnts]
    area = pcb.AddArea(None, net, pcb.GetLayerID(layer), pnt.to_VEC2I(pnts[0]), pcbnew.ZONE_BORDER_DISPLAY_STYLE_DIAGONAL_EDGE)
    poly = area.Outline()
    for idx, pt in enumerate(pnts):
        if idx == 0:
            continue
        poly.Append(pt[0], pt[1])
    return area


def add_zone(rect, layer, net_name="GND"):
    area = _add_area(rect, layer, net_name)
    return area


def add_rule_area(pnts, layer):
    area = _add_area(pnts, layer, None)
    area.SetIsRuleArea(True)
    area.SetDoNotAllowTracks(False)
    area.SetDoNotAllowVias(False)
    area.SetDoNotAllowPads(False)
    area.SetDoNotAllowCopperPour(True)


def make_rect(size, offset=(0, 0)):
    FourCorners = [(0, 0), (1, 0), (1, 1), (0, 1)]
    return [vec2.add((size[0] * cnr[0], size[1] * cnr[1]), offset) for cnr in FourCorners]
