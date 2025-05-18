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
                    ("C1", (+8.0, 5.66 - 0.89 * 2), 0),
                    ("C2", (+8.0, 5.66 - 0.89 * 4), 0),
                    ("C5", (+8.0, 5.66 - 0.89 * 6), 0),
                    ("C3", (+8.0, 5.66 - 0.89 * 8), 0),
                    ("R1", (-8.0, 5.66 - 0.89 * 3), 0),
                    ("C6", (-8.0, 5.66 - 0.89 * 1), 0),
                    ("R2", (-8.0, 5.66 - 0.89 * 9), 0),
                    ("C4", (-8.0, 5.66 - 0.89 * 7), 0),
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
                (13.2, 9),
                90,
                [
                    ("U2", (0, 0), 0),
                    ("C11", (-2.8, 0), -90),
                    ("C12", (-4.3, 0), -90),
                    ("C13", (+2.8, 0), -90),
                    ("R3", (+4.3, 0), -90),
                    ("R4", (+0.9, +2.5), 180),
                    ("C14", (+0.9, -2.5), 180),
                ],
            ),
        ],
    )
    for n in range(5):
        kad.move_mods(
            vec2.add(board_orig, (-13, -2.54 * 2.5 * (n - 2))),
            90,
            [
                (f"R{11+2*n}", (-1.5, 0.8), 90),
                (f"R{12+2*n}", (0, 0.8), -90),
                (f"C2{n+1}", (+1.5, 0.8), 90),
                (f"J{n+1}", (-1.27, -2.54), 90),
            ],
        )
    for iy, dy in enumerate((-1, +1)):
        for ix, dx in enumerate((-1, +1)):
            idx = 2 * iy + ix + 1
            kad.set_mod_pos_angle(f"H{idx}", vec2.add(board_orig, ((board_width / 2 - 3) * dx, (board_height / 2 - 3) * dy)), 0)


w_pwr, r_pwr = 0.50, 1.5  # power
w_dat, r_dat = 0.30, 0.7  # row / col

r_tri = 0.7


def wire_mod():
    pmw = "U1"
    xiao = "U3"
    ldo = "U2"
    kad.wire_mod_pads(
        [
            ### PMW3360
            # ("C1", "2", "C3", "2", w_dat, (Strt), "F.Cu"),
            ("C1", "1", "C2", "1", w_dat, (Strt)),
            (pmw, "3", "C2", "1", w_dat, (Strt)),
            (pmw, "4", "C5", "1", w_dat, (Strt)),
            (pmw, "5", "C3", "1", w_dat, (Strt)),
            (pmw, "12", "R2", "2", w_dat, (Strt)),
            (pmw, "13", "C4", "2", w_dat, (Strt)),
            ("R2", "1", "C4", "1", w_dat, (Dird, 0, 90)),
            (pmw, "15", "R1", "2", w_dat, (Strt)),
            ("C6", "1", "R1", "1", w_dat, (Dird, 0, 90)),
            ### LDO
            (ldo, "1", "C11", "1", w_dat, (Dird, 0, -45)),
            (ldo, "2", "C11", "2", w_dat, (Dird, 0, -45)),
            ("C11", "1", "C12", "1", w_dat, (Strt)),
            ("C11", "2", "C12", "2", w_dat, (Strt)),
            ("C14", "2", "R4", "2", w_dat, (Dird, [(0, 0.1), 90], 0)),
            ("C14", "2", ldo, "2", w_dat, (Dird, [(0, 0.1), 90], 0, 0)),
            # ("C14", "2", ldo, "2", w_dat, (Dird, [(0, 0.1), (90, 8), 90], 0, r_dat)),
        ]
    )

    return
    via_3v3_led = kad.add_via(kad.calc_pos_from_pad(rj45, "23", (0.6, -1.8)), VCC, via_size_pwr)
    # RJ45 - xiao_r
    via_5vd = kad.add_via_relative(rj45, "2", (0, -3.8), via_size_pwr)
    kad.wire_mod_pads(
        [
            # 5VD
            (xiao_r, "1", rj45, via_5vd, w_pwr, (Dird, 45, 90), "F.Cu"),
            # 5VD
            (rj45, "4", rj45, via_5vd, w_pwr, (Dird, [(+90, 1.5), 0], 90), "F.Cu"),
            (rj45, "16", rj45, via_5vd, w_pwr, (Dird, [(+90, 1.5), 0], 90, r_led), "F.Cu"),
            (rj45, "16", rj45, "4", w_pwr, (Dird, [(+90, 1.5), 0], 90, r_led), "F.Cu"),
            # GND
            # (xiao_r, "2", rj45, "2", w_pwr, (Dird, -45, 90), "In2.Cu"),
            # 3V3
            (xiao_r, "3", rj45, "6", w_pwr, (ZgZg, 0, 30), "In2.Cu"),
            # 3V3
            (rj45, "6", rj45, "18", w_pwr, (Dird, [(+90, 1.5), 0], 90, r_tri), "B.Cu"),
            (rj45, "18", rj45, via_3v3_led, w_pwr, (Dird, 90, -45, r_pwr), "F.Cu"),
            # LED1/2/3
            (xiao_r, "4", rj45, "9", w_led, (Dird, -45, [(-135, 2.0), (180, 1.7), (135, 2.0), 0], r_led), "In1.Cu"),
            (xiao_r, "5", rj45, "11", w_led, (Dird, -45, 90), "In1.Cu"),
            (xiao_r, "7", rj45, "21", w_led, (Dird, -45, 90), "In1.Cu"),
        ]
    )
    # RJ45 - xiao_l
    via_sda = kad.add_via_relative(rj45, "19", (2.0, 2.4), via_size_dat)
    via_sck = kad.add_via_relative(rj45, "15", (2.0, 2.4), via_size_dat)
    kad.wire_mod_pads(
        [
            # Full colr LED
            (xiao_l, "1", rj45, "1", w_dat, (Dird, [(0, 2), -45], 90), "In1.Cu"),
            (xiao_l, "2", rj45, "13", w_dat, (Dird, [(0, 3.2), 90], [(0, 5.6), -45], r_led), "In1.Cu"),
            # NRST
            (xiao_l, "3", rj45, "5", w_dat, (ZgZg, 0, 45), "F.Cu"),
            (xiao_l, "3", rj45, "17", w_dat, (Dird, [(0, 2.2), 90], 90, r_led), "In2.Cu"),
            # SCK
            (rj45, "3", rj45, "15", w_dat, (Dird, [(-90, 2.4), 0], 90, r_dat), "B.Cu"),
            (rj45, "3", xiao_l, "6", w_dat, (Dird, [(-90, 2.4), 0], 0, r_dat), "B.Cu"),
            (rj45, "15", xiao_l, "6", w_dat, (Dird, [(-90, 2.4), 0], 0, r_dat), "B.Cu"),
            (rj45, "3", "R12", "1", w_dat, (Dird, [(-90, 2.4), 0], 90, r_dat), "B.Cu"),
            (rj45, "15", "R12", "1", w_dat, (Strt), "B.Cu"),
            # SDA
            (rj45, "7", rj45, "19", w_dat, (Dird, [(-90, 2.4), 0], 90, r_dat), "F.Cu"),
            (rj45, "7", xiao_l, "5", w_dat, (Dird, [(-90, 2.4), 0], 0, r_dat), "F.Cu"),
            (rj45, "19", xiao_l, "5", w_dat, (Dird, [(-90, 2.4), 0], 0, r_dat), "F.Cu"),
            (rj45, "19", "R11", "1", w_dat, (Dird, 90, 0, r_dat), "B.Cu"),
        ]
    )
    # I2C 3V3
    via_left = kad.add_via_relative(rj45, "20", (-1.6, 0), via_size_pwr)
    kad.wire_mod_pads(
        [
            (rj45, "18", rj45, via_left, w_led, (Dird, [(90, 1.5), 0], 90, r_led), "B.Cu"),
            (rj45, via_left, "R12", "2", w_led, (Dird, [(-90, 0.4), (-60, 3.2), 90], 0, r_led), "B.Cu"),
            ("R11", "2", "R12", "2", w_led, (ZgZg, 90, 45), "B.Cu"),
        ]
    )
    # LED
    kad.wire_mod_pads(
        [
            # LED1/2/3
            (rj45, "24", "R4", "2", w_led, (Dird, 90, [(0, 1.6), 90], r_led), "B.Cu"),
            (rj45, "22", "R3", "2", w_led, (Dird, 90, [(0, 0.0), 90], r_led), "B.Cu"),
            (rj45, "12", "R2", "2", w_led, (Dird, 90, [(0, 1.0), 90], r_led), "B.Cu"),
            (rj45, "10", "R1", "2", w_led, (Dird, 0, [(0, 1.8), (-90, 10), (-135, 0.5 * 1.414), (-90, 3.2), -135], r_led), "B.Cu"),
            # LED4
            (rj45, "23", rj45, via_left, w_led, (ZgZg, 90, 30), "In1.Cu"),
            (xiao_l, "4", rj45, via_left, w_led, (Dird, [(0, 2.2), 90], [(-90, 0.4), -60], r_led), "In1.Cu"),
            # 3V3
            (rj45, via_3v3_led, "R4", "1", w_pwr, (Dird, -45, 90, r_tri), "B.Cu"),
            ("R4", "1", "R3", "1", w_pwr, (Strt), "B.Cu"),
            ("R3", "1", "R4", "1", w_pwr, (Strt), "B.Cu"),
        ]
    )
    # C1/2
    kad.wire_mod_pads(
        [
            ("C1", "2", rj45, "2", w_pwr, (Dird, 45, 90), "B.Cu"),
            ("C1", "1", rj45, "4", w_pwr, (Dird, 135, 90), "B.Cu"),
            ("C2", "2", rj45, "8", w_pwr, (Dird, 45, 0), "B.Cu"),
        ]
    )
    for via in [via_5vd, via_left, via_sda, via_sck]:
        pcb.Delete(via)

    # gnd vias
    for pad in "1234567":
        kad.add_via(kad.calc_pos_from_pad(xiao_l, pad, (1.2, -2.54 / 2)), GND, via_size_dat)
    for pad in ["1", "3", "5"]:
        kad.add_via(kad.calc_pos_from_pad(rj45, pad, (-2.02 / 2, 1.4)), GND, via_size_dat)
    for pad in ["15", "17"]:
        kad.add_via(kad.calc_pos_from_pad(rj45, pad, (+2.02 / 2, 1.4)), GND, via_size_dat)
    kad.add_via(kad.calc_pos_from_pad(rj45, "24", (+2.29 / 2, +1.2)), GND, via_size_dat)
    kad.add_via(kad.calc_pos_from_pad(rj45, "24", (+2.29 / 2, -1.2)), GND, via_size_dat)
    kad.add_via(kad.calc_pos_from_pad(rj45, "10", (+2.29 / 2, -1.2)), GND, via_size_dat)
    kad.add_via(kad.calc_pos_from_pad(rj45, "20", (-2.8, 0)), GND, via_size_dat)
    kad.add_via(kad.calc_pos_from_pad(rj45, "20", (0, 6.4)), GND, via_size_dat)
    kad.add_via(kad.calc_pos_from_pad(rj45, "11", (1.8, 0)), GND, via_size_dat)


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
        (2.2, 180 - 18, -90, [f"R{2*n+11}" for n in range(5)]),
        (2.2, 18, -90, [f"R{2*n+12}" for n in range(5)]),
        (1.3, -90, 0, [f"C2{n+1}" for n in range(5)]),
        (0.7 + 1.5 * 2, -90, 0, ["C25"]),
        (1.0 + 1.5 * 2, +90, 0, ["R20"]),
        (1.3 + 1.5 * 2, -90, 0, ["R19"]),
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
        (2, "NC"),
        (3, "VPIX"),
        (4, "1V9"),
        (5, "3V3"),
        (6, "NC"),
        (7, "nRST"),
        (8, "GND"),
        (9, "MOTN"),
        (10, "SCLK"),
        (11, "MOSI"),
        (12, "MISO"),
        (13, "GND"),
        (14, "NC"),
        (15, "LED"),
        (16, "NC"),
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
