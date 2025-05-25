# exec(open('place_and_route.py').read())
import importlib

import pcbnew
from kadpy import kad, mat2, pnt, vec2

importlib.reload(kad)
importlib.reload(pnt)
importlib.reload(vec2)
importlib.reload(mat2)

kad.UnitMM = True
kad.PointDigits = 3

# alias
Strt = kad.Straight
Dird = kad.Directed
ZgZg = kad.ZigZag

Round = kad.Round
BezierRound = kad.BezierRound
LinearRound = kad.LinearRound

# in mm
VIA_Size = [(1.2, 0.6), (1.15, 0.5), (0.92, 0.4), (0.8, 0.3)]

via_size_pwr = VIA_Size[1]
via_size_dat = VIA_Size[2]
via_size_gnd = VIA_Size[3]

Cu_layers = ["F.Cu", "B.Cu"]  # , "In1.Cu", "In2.Cu"]

pcb = pcbnew.GetBoard()
# for n in dir(pcb):
#     if "Zone" in n:
#         print(n)

GND = pcb.FindNet("GND")
VCC = pcb.FindNet("3V3")

board_width = 34
board_height = 44  # 2.54 * 7
board_size = (board_width, board_height)
board_orig = (100, 100)


def place_mods():
    kad.move_mods(
        board_orig,
        0,
        [
            (  # pmw3360
                None,
                (0, -10.5),
                180,
                [
                    ("U1", (0, 0), 0),
                    ("C1", (+8.0, 5.66 - 0.89 * 1), 0),
                    ("C2", (+8.0, 5.66 - 0.89 * 3), 0),
                    ("C5", (+8.0, 5.66 - 0.89 * 5), 0),
                    ("C6", (+8.0, 5.66 - 0.89 * 7), 0),
                    ("C3", (+8.0, 5.66 - 0.89 * 9), 0),
                    ("C4", (+8.0, 5.66 - 0.89 * 11), 0),
                    ("C7", (-8.0, 5.66 - 0.89 * -1), 0),
                    ("C8", (-8.0, 5.66 - 0.89 * 1), 0),
                    ("R1", (-8.0, 5.66 - 0.89 * 3), 0),
                    ("JP1", (-8.0, 5.66 - 0.89 * 6), 180),
                    ("R2", (-8.0, 5.66 - 0.89 * 9), 0),
                    ("C9", (-8.0, 5.66 - 0.89 * 11), 0),
                ],
            ),
            (  # xiao
                None,
                (0, +11),
                180,
                [
                    ("U3", (-0.20, 0), 180),
                    ("J7", (+2.54 * 3, -2.54 * 3), 0),
                    ("J8", (-2.54 * 3, -2.54 * 3), 0),
                    ("J9", (0.54, 5.6), 180),
                ],
            ),
            (  # LDO
                None,
                (13, 10),
                90,
                [
                    ("U2", (0, 0), 0),
                    ("C11", (-2.8, 0), -90),
                    ("C12", (-4.3, 0), -90),
                    ("C13", (+5.8, 0), -90),
                    ("R3", (+4.3, 0), -90),
                    ("R4", (+2.8, 0), +90),
                    ("C14", (+0.0, -2.5), 180),
                ],
            ),
        ],
    )
    # SW
    for n in range(5):
        kad.move_mods(
            vec2.add(board_orig, (-2.54 * 5, -2.54 * 2.4 * (n - 2) + 1.0)),
            90,
            [
                (f"R{11+2*n}", (+1.5, 0.74), 90),
                (f"R{12+2*n}", (0, 0.72), -90),
                (f"C2{n+1}", (-1.5, 0.72), 90),
                (f"J{n+1}", (+1.27, -2.54), -90),
            ],
        )
    for iy, dy in enumerate((-1, +1)):
        for ix, dx in enumerate((-1, +1)):
            idx = 2 * iy + ix + 1
            kad.set_mod_pos_angle(f"H{idx}", vec2.add(board_orig, ((board_width / 2 - 3) * dx, (board_height / 2 - 3) * dy)), 0)


w_pwr, r_pwr = 0.40, 0.8  # power
w_dat, r_dat = 0.32, 0.6  # data

r_tri = 0.7


def wire_mod():
    pmw = "U1"
    ldo = "U2"
    xiao = "U3"

    via_ldo_vbus = [
        kad.add_via_relative("C12", "1", (-1.4, 0), via_size_pwr),
        kad.add_via_relative(ldo, "3", (-1.2, +1.2), via_size_dat),
    ]

    via_gnd_U1 = kad.add_via(vec2.add(kad.get_mod_pos(pmw), (0, -10)), GND, via_size_pwr)
    via_gnd_C1 = kad.add_via_relative("C1", "2", (0, +1.6), via_size_pwr)
    # via_gnd_C4 = kad.add_via_relative("C4", "2", (0, +1.6), via_size_pwr)
    # via_gnd_C6 = kad.add_via_relative("C6", "2", (0, +1.6), via_size_pwr)
    via_1V9_R1 = kad.add_via_relative("R1", "1", (-1.6, 0), via_size_pwr)
    via_1V9_R3 = kad.add_via_relative("R3", "1", (-1.0, 0), via_size_pwr)

    kad.wire_mod_pads(
        [
            ### PMW3360
            # left
            (pmw, "3", "C2", "1", w_dat, (Dird, -45, 0, r_dat)),  # VDDPIX
            ("C1", "1", "C2", "1", w_dat, (Strt)),
            (pmw, "4", "C5", "1", w_dat, (Dird, 0, 90, 0)),  # 1V9
            (pmw, "4", "C6", "1", w_dat, (Dird, 0, 90, 0)),  # 1V9
            (pmw, "5", "C3", "1", w_dat, (Dird, +45, 0, r_dat)),  # 3V4
            ("C3", "1", "C4", "1", w_dat, (Strt)),
            (pmw, "8", "C4", "2", w_pwr, (Dird, [(0, 1.2), 90], [(90, 1.2), 0], r_pwr)),  # GND
            ("C1", "2", "C4", "2", w_pwr, (Strt)),
            # right
            ("C7", "1", "C8", "1", w_dat, (Strt)),
            ("C7", "2", "C8", "2", w_dat, (Strt)),
            ("R1", "1", None, via_1V9_R1, w_dat, (Strt)),
            ("R1", "1", "C8", "1", w_dat, (Dird, 90, 0)),
            (pmw, "15", "R1", "2", w_dat, (Strt)),  # LED_P
            (pmw, "13", "JP1", "1", w_dat, (Dird, 45, 0, r_dat)),  # nCS
            ("R2", "1", "C9", "1", w_dat, (Dird, 0, 90)),
            (pmw, "12", "R2", "2", w_dat, (Strt)),  # MISO
            # (pmw, "13", "C4", "2", w_dat, (Strt)),
            # ("C1", via_gnd_C1, "C3", via_gnd_C3, w_dat, (Strt), "F.Cu"),
            # ("C1", via_gnd_C1, "C3", via_gnd_C3, w_dat, (Strt), "F.Cu"),
            ("C1", via_gnd_C1, None, via_gnd_U1, w_dat, (Dird, 90, 0, r_pwr)),
            # ("C6", via_gnd_C6, None, via_gnd_U1, w_dat, (Dird, 90, 0, r_pwr)),
            # ("C4", "2", None, via_gnd_C4, w_dat, (Strt)),
            # ("C6", "2", None, via_gnd_C6, w_dat, (Strt)),
            ### LDO
            ("C12", "1", None, via_ldo_vbus[0], w_dat, (Strt)),
            (ldo, "3", None, via_ldo_vbus[1], w_dat, (Dird, 90, 0, 0.4)),
            (ldo, via_ldo_vbus[0], None, via_ldo_vbus[1], w_dat, (Dird, 90, 0, r_dat), "B.Cu"),
            (ldo, "1", "C11", "1", w_dat, (Dird, 0, -45)),
            (ldo, "2", "C11", "2", w_dat, (Dird, 0, -45, 0.4)),
            (ldo, "2", "R4", "2", w_dat, (Dird, 0, -45, 0.4)),
            (ldo, "2", "C14", "2", w_dat, (Dird, 0, [(135, 1.1), 90], 0)),
            (ldo, "4", "R4", "1", w_dat, (Dird, 0, 0, 0)),
            (ldo, "5", "C14", "1", w_dat, (Dird, 0, 90, 0)),
            ("R3", "2", "R4", "1", w_dat, (Strt)),
            ("R3", "1", "C13", "1", w_dat, (Dird, 90, 0)),
            ("R3", "2", "C13", "2", w_dat, (Dird, 90, 0)),
            ("C11", "1", "C12", "1", w_dat, (Strt)),
            ("C11", "2", "C12", "2", w_dat, (Strt)),
            ("R3", "1", None, via_1V9_R3, w_dat, (Strt)),
            ("C14", "1", "R3", via_1V9_R3, w_dat, (Dird, [(180, 1.0), 135], 90, 0.4)),
            ### LDO <--> pmw
            ("R1", via_1V9_R1, "R3", via_1V9_R3, w_pwr, (Dird, 0, 90, 0.4)),
            ("R1", via_1V9_R1, "U1", "4", w_pwr, (Dird, [(-90, 3.4), (-45, 4.0), 0], [(0, 1.6), 90], r_pwr), "In2.Cu"),
            ### LD0 <--> xiao
            (xiao, "13", "C12", "2", w_pwr, (Dird, 0, 90, r_pwr), "F.Cu"),  # GND
            (xiao, "14", "C12", via_ldo_vbus[0], w_pwr, (Dird, 0, 90, r_pwr), "B.Cu"),  # VBUS
            ### xiao <--> pmw
            (xiao, "7", pmw, "7", w_dat, (Dird, 90, -135, r_dat)),  # nRESET
            (xiao, "8", pmw, "9", w_dat, (Dird, 90, -45, r_dat), "B.Cu"),  # Motion
            (xiao, "9", pmw, "10", w_dat, (Dird, [(0, 2.4), 90], -45, r_dat), "B.Cu"),  # SCLK
            (xiao, "10", pmw, "12", w_dat, (Dird, [(0, 3.4), 90], -45, r_dat), "In2.Cu"),  # MISO
            (xiao, "11", pmw, "11", w_dat, (Dird, [(0, 3.4), 90], -45, r_dat), "B.Cu"),  # MOSI
            ("J8", "3", "C3", "1", w_pwr, (Dird, [(-30, 2.8), 90], [(90, 9), 0], r_pwr), "F.Cu"),  # 3V3
            ("J8", "3", "R2", "1", w_pwr, (Dird, [(-30, 2.8), 90], [(90, 8.11), 0], r_pwr), "F.Cu"),  # 3V3
            ("C3", "1", "R2", "1", w_pwr, (Dird, 90, [(90, 8.11), 0], r_pwr), "F.Cu"),
        ]
    )
    pcb.Remove(via_gnd_U1)
    pcb.Remove(via_1V9_R3)

    # SW
    via_sw_3v3 = [kad.add_via_relative(f"R{11+2*n}", "1", (0.2, 1.5), via_size_pwr) for n in range(5)]
    for n in range(5):
        kad.wire_mod_pads(
            [
                (f"J{n+1}", "1", f"R{11+2*n}", "2", w_dat, (Dird, 90, 0, 0)),
                (f"J{n+1}", "2", f"C{21+n}", "2", w_dat, (Dird, 90, 0, 0)),
                (f"R{12+2*n}", "2", f"C{21+n}", "1", w_dat, (Dird, 90, 0, 0)),
                (f"R{11+2*n}", "2", f"R{12+2*n}", "1", w_dat, (Dird, 0, 90, 0)),
                (f"R{11+2*n}", "1", None, via_sw_3v3[n], w_dat, (Dird, 90, 0, 0)),
            ]
        )

    # SW --> xiao
    r_sw = 1.2
    kad.wire_mod_pads(
        [
            ("C21", "1", "J7", "2", w_dat, (Dird, 0, -45, r_sw)),
            ("C22", "1", "J7", "3", w_dat, (Dird, [(180, 1.1), 90], -45, r_sw)),
            ("C23", "1", "J7", "4", w_dat, (Dird, [(180, 1.1), (90, 5.2), (135, 1.2), 90], -45, r_sw)),
            ("C24", "1", "J7", "5", w_dat, (Dird, [(180, 1.1), (90, 5.2), (135, 1.2), 90], -45, r_sw)),
            ("C25", "1", "J7", "6", w_dat, (Dird, [(180, 1.1), (90, 5.2), (135, 1.2), (90, 4.6), (135, 1.2), 0], [(-135, 2.2), -90], r_sw)),
        ]
    )

    # SW 3V3 pwr
    via_sw_vert = kad.add_via_relative("R11", "1", (-1.0, -4.4 - 1.0), via_size_pwr)
    via_xiao_btm = kad.add_via_relative(xiao, "13", (-8.0, -2.54 * 1.8), via_size_pwr)
    kad.wire_mod_pads(
        [
            ("J1", via_sw_vert, "J5", via_sw_3v3[4], w_pwr, (Dird, -45, 90, r_pwr), "In1.Cu"),  # vertical 3V3
            ("J1", via_sw_vert, None, via_xiao_btm, w_pwr, (Dird, [(-45, 1.2), 90], 0, r_pwr), "In1.Cu"),  # 3V3
            ("J1", via_sw_vert, None, via_xiao_btm, w_pwr, (Dird, [(-45, 1.2), 90], 0, r_pwr), "F.Cu"),  # GND
            ("J8", "3", None, via_xiao_btm, w_pwr, (Dird, [(30, 2.8), 90], 0, r_pwr), "In1.Cu"),  # 3V3
            ("J8", "2", None, via_xiao_btm, w_pwr, (Dird, [(30, 2.8), 90], 0, r_pwr), "F.Cu"),  # GND
            ("J1", "2", None, via_sw_vert, w_pwr, (Dird, [(-90, 1.6), 0], -45, r_pwr), "F.Cu"),  # GND
        ]
    )
    pcb.Remove(via_sw_vert)
    pcb.Remove(via_xiao_btm)


def draw_edge_cuts():
    width = 0.12

    Radius = 2.6
    _org = board_orig
    cnrs = [
        ((vec2.add(_org, (0, -board_height / 2)), 0), Round, [Radius]),
        ((vec2.add(_org, (+board_width / 2, 0)), 90), Round, [Radius]),
        ((vec2.add(_org, (0, +board_height / 2)), 180), Round, [Radius]),
        ((vec2.add(_org, (-board_width / 2, 0)), 270), Round, [Radius]),
    ]
    kad.draw_closed_corners(cnrs, "Edge.Cuts", width)


# References
def set_text_prop(text, pos, angle, offset_length, offset_angle, text_angle):
    if text_angle == None:
        text.SetVisible(False)
    else:
        text.SetVisible(True)
        # tsz = 1.0
        tsz = 0.9
        text.SetTextSize(pnt.to_VEC2I(pcbnew.wxSizeMM(tsz, tsz)))
        text.SetTextThickness(pcbnew.FromMM(0.18))
        pos_text = vec2.scale(offset_length, vec2.rotate(-(offset_angle + angle)), pos)
        text.SetPosition(pnt.to_VEC2I(pnt.to_unit(vec2.round(pos_text, 3), True)))
        text.SetTextAngle(kad.to_ANGLE(text_angle))
        text.SetKeepUpright(False)


def set_refs():
    # hide value texts
    for mod in pcb.GetFootprints():
        ref = mod.Reference()
        val = mod.Value()
        val.SetVisible(False)
    # hide mounting holes
    for n in range(4):
        mod = kad.get_mod(f"H{n+1}")
        ref = mod.Reference()
        ref.SetVisible(False)
    refs = [
        # PMW3360
        (9.6, -90, 0, ["U1"]),
        (2.0, 180, -90, ["R1", "R2", "C4", "C6"]),
        (2.0, 0, +90, ["C1", "C2", "C3", "C5"]),
        # LDO
        (0, 0, 0, ["U2"]),
        (3.0, +90, 0, ["R3", "C13"]),
        (3.0, -90, 0, ["C11", "C12"]),
        (3.2, 180, 90, ["R4", "C14"]),
        # xiao
        (14, 55, 0, ["U3"]),
        (None, None, None, ["J7", "J8", "J9"]),
        # switches
        (1.8, 90, 0, [f"J{n+1}" for n in range(5)]),
        (2.2, 180 + 18, -90, [f"R{2*n+11}" for n in range(5)]),
        (2.2, -18, -90, [f"R{2*n+12}" for n in range(5)]),
        (1.3, 90, 0, [f"C2{n+1}" for n in range(5)]),
        (1.3, -90, 0, ["R17"]),
        (2.0, 0, -90, ["R18"]),
        (1.3 + 1.5 * 2, -90, 0, ["C25"]),
        (1.0 + 1.5 * 2, +90, 0, ["R20"]),
        (0.7 + 1.5 * 2, -90, 0, ["R19"]),
    ]
    for offset_length, offset_angle, text_angle, mod_names in refs:
        for mod_name in mod_names:
            mod = kad.get_mod(mod_name)
            if mod is None:
                continue
            pos, angle = kad.get_mod_pos_angle(mod_name)
            ref = mod.Reference()
            set_text_prop(ref, pos, angle, offset_length, offset_angle, text_angle)
            # for item in mod.GraphicalItems():
            #     print(f"{type(item) = }")
            # if type(item) is pcbnew.PCB_TEXT and item.GetShownText() == ref.GetShownText():
            #     set_text_prop(item, pos, angle, offset_length, offset_angle, text_angle)
    tsz = 0.9
    # J7
    angle = 180
    pads = ["DSW1", "DSW2", "DSW3", "DSW4", "DSW5", "nCS", "nRST"]
    for idx, pad in enumerate(pads):
        pos = kad.calc_pos_from_pad("J7", f"{idx+1}", (-1.6, 0))
        kad.add_text(pos, angle, pad, "F.SilkS", (tsz, tsz), 0.18, pcbnew.GR_TEXT_H_ALIGN_RIGHT, pcbnew.GR_TEXT_V_ALIGN_CENTER)
        kad.add_text(pos, angle, pad, "B.SilkS", (tsz, tsz), 0.18, pcbnew.GR_TEXT_H_ALIGN_LEFT, pcbnew.GR_TEXT_V_ALIGN_CENTER)
    # J8
    angle = 180
    pads = ["VBUS", "GND", "3V3", "MOSI", "MISO", "SCLK", "MOTN"]
    for idx, pad in enumerate(pads):
        pos = kad.calc_pos_from_pad("J8", f"{idx+1}", (+1.6, 0))
        kad.add_text(pos, angle, pad, "F.SilkS", (tsz, tsz), 0.18, pcbnew.GR_TEXT_H_ALIGN_LEFT, pcbnew.GR_TEXT_V_ALIGN_CENTER)
        kad.add_text(pos, angle, pad, "B.SilkS", (tsz, tsz), 0.18, pcbnew.GR_TEXT_H_ALIGN_RIGHT, pcbnew.GR_TEXT_V_ALIGN_CENTER)
    # J9
    angle = 90
    pads = ["BAT+", "BAT-"]
    for idx, pad in enumerate(pads):
        pos = kad.calc_pos_from_pad("J9", f"{idx+1}", (0, -1.6))
        kad.add_text(pos, angle, pad, "F.SilkS", (tsz, tsz), 0.18, pcbnew.GR_TEXT_H_ALIGN_LEFT, pcbnew.GR_TEXT_V_ALIGN_CENTER)
        kad.add_text(pos, angle, pad, "B.SilkS", (tsz, tsz), 0.18, pcbnew.GR_TEXT_H_ALIGN_RIGHT, pcbnew.GR_TEXT_V_ALIGN_CENTER)
    # U1
    angle = kad.get_mod_angle("U1")
    pads = [
        (2, "nc"),
        (3, "VPIX"),
        (4, "1V9"),
        (5, "3V3"),
        (6, "nc"),
        (7, "nRST"),
        (8, "GND"),
        (9, "MOTN"),
        (10, "SCLK"),
        (11, "MOSI"),
        (12, "MISO"),
        (13, "GND"),
        (14, "nc"),
        (15, "LED"),
        (16, "nc"),
    ]
    for idx, (pin, pad) in enumerate(pads):
        side = +1 if pin < 9 else -1
        pos = kad.calc_pos_from_pad("U1", f"{pin}", (side, 0))
        if True:  # Front
            _pos = pos
            if pin in [2, 3, 4, 5]:
                _pos = vec2.add(_pos, (0, -1.78 * 4 - (pin - 5) * 0.4))
            elif pin in [12, 13, 15, 16]:
                _pos = vec2.add(_pos, (4.4, 0))
            halign = pcbnew.GR_TEXT_H_ALIGN_LEFT if side > 0 else pcbnew.GR_TEXT_H_ALIGN_RIGHT
            kad.add_text(_pos, angle, pad, "F.SilkS", (tsz, tsz), 0.18, halign, pcbnew.GR_TEXT_V_ALIGN_CENTER)
        halign = pcbnew.GR_TEXT_H_ALIGN_LEFT if side < 0 else pcbnew.GR_TEXT_H_ALIGN_RIGHT
        kad.add_text(pos, angle, pad, "B.SilkS", (tsz, tsz), 0.18, halign, pcbnew.GR_TEXT_V_ALIGN_CENTER)


def add_zone(net_name, layer_name, rect):
    settings = pcb.GetZoneSettings()
    settings.m_ZoneClearance = pcbnew.FromMils(12)
    pcb.SetZoneSettings(settings)

    zone = kad.add_zone(rect, layer_name, net_name)
    zone.SetMinThickness(pcbnew.FromMils(13))
    zone.SetThermalReliefGap(pcbnew.FromMils(13))
    # zone.Hatch()


def main():
    place_mods()
    wire_mod()
    draw_edge_cuts()
    set_refs()

    # logo
    for mod, angle in [("G1", 180)]:
        if kad.get_mod(mod) is not None:
            kad.move_mods((90, 87), 0, [(mod, (0, 0), angle)])

    # zones
    rect = kad.make_rect(vec2.scale(1.1, board_size), vec2.scale(-0.5 * 1.1, board_size, board_orig))
    # for layer in Cu_layers:
    #     add_zone("GND", layer, rect)

    # name
    kad.add_text(
        (board_orig[0] + 16, board_orig[1] - 6),
        90,
        f"orihikarna 2025/06/30",
        "F.Silkscreen",
        (0.8, 0.8),
        0.4,
        pcbnew.GR_TEXT_H_ALIGN_CENTER,
        pcbnew.GR_TEXT_V_ALIGN_CENTER,
    )


if __name__ == "__main__":
    kad.removeDrawings()
    kad.removeTracksAndVias()
    main()
    pcbnew.Refresh()
