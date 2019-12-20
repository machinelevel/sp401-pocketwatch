"""
This file is Copyright 2019 Eric R. Johnston
Patent pending

LICENSE:
This software may be used only with written permission from the author.
If you want to use it, please email me at pocketwatch@machinelevel.com

"""
import numpy as np
from ej_outer_tooth_shapes import outer_tooth_shapes_p10
try:
    import collada as col  # http://pycollada.github.io/creating.html
    have_collada = True
except:
    have_collada = False

# Global scale values
# All meassurements are in mm
use_collada = False
strut_outer_radius = 1.0
thinnest_material_wall = 0.4
spur_teeth_thickness = thinnest_material_wall
min_material_thickness = 0.4
slide_buffer_dist = 0.1
spur_tooth_pitch = 2.0
tooth_rim_thickness = 0.4
geneva_thickness = 1.0
each_platter_z = tooth_rim_thickness * 2 * 4.0
geneva_outer_radius = 20.0
planet_outer_ref_radius = 15.0
ball_bearing_radius = 0.5 * thinnest_material_wall


strut_outer_radius = 1.0
strut_inner_radius = strut_outer_radius - thinnest_material_wall

"""
      best_errors [0.0002094545855868546, 0.002124126164289919, 0.0020335698811777547, 0.00014134032982582312]
      best_ratios [231.75420417869884, 1167.8433179723504, 1559.8970466919648, 59.061319340329824]
      best_gears M[(-27, 29, 13, 29, 23, -13, 21, 39),
                 V (-27, 29, 13, 31, 38, -29, 14, 38), 
                 M (-27, 29, 13, 11, 26, -28, 19, 37), 
                 L (-27, 29, 13, 23, 13, -7, 30, 37)]
      error_degrees_per_year [0.2376797302003289, 0.4783273624581241, -0.34284199857996345, 0.6293529646316143]
"""

ball_bearing_faces = None

all_platters = {
    'Drive': {
        'use_dual_drive': False,
        'pos': [0.0, 0.0, -2 * each_platter_z],
        'planetary_mult': 1,
        'scale_geneva': 1.0,
        'scale_planetary': 1.73 * 0.95,
        'rings_under_planets':False,
        'support_platter':None,
        'rotor_azimuth':-90,
        'scale_base_ring_radius':1.0,
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':27},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':13},
            'Pr'     :{'type':'spur',   'inout':'in', 'teeth':29, 'outer_rail':3.5},
            'Grotor' :{'type':'rotor'},
            'feet'   :{'type':'feet'},
            'shaft'  :{'type':'shaft'},
        },
    },
    'Mars': {
        'use_dual_drive': False,
        'pos': [0.0, 0.0, -1 * each_platter_z],
        'planetary_mult': 1,
        'scale_geneva': 1.0,
        'scale_planetary': 2.4,
        'rings_under_planets':True,
        'support_platter':None,
        'rotor_azimuth':-90,
        'scale_base_ring_radius':1.0,
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':28},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':19, 'inner_rail':0},
            'Pr'     :{'type':'spur',   'inout':'in', 'teeth':37, 'outer_rail':3.5},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':11},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':26, 'outer_rail':1},
            'Dr1'     :{'type':'spur',   'inout':'out', 'teeth':int(26 * 0.33)},
            'Dr2'     :{'type':'spur',   'inout':'out', 'teeth':int(26 * 0.33)},
            'Grotor' :{'type':'rotor'},
            'feet'   :{'type':'feet'},
            'shaft'  :{'type':'shaft'},
        },
    },
    'Luna': {
        'use_dual_drive': False,
        'pos': [0.0, 0.0, 0 * each_platter_z],
        'planetary_mult': 1,
        'scale_geneva': 0.75,
        'scale_planetary': 1.0,
        'support_platter':None,
        'rotor_azimuth':180,
        'scale_base_ring_radius':1.0,
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':7*2},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':30, 'inner_rail':0},
            'Pp0'    :{'type':'spur',   'inout':'out', 'teeth':37, 'p_placement':(0.0 * np.pi, 0), 'roll_to':(90, 'Ps')},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':23*2},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':13, 'outer_rail':1},
            'Dr1'     :{'type':'spur',   'inout':'out', 'teeth':int(13 * 0.33)},
            'Dr2'     :{'type':'spur',   'inout':'out', 'teeth':int(13 * 0.33)},
            'Grotor' :{'type':'rotor'},
            'feet'   :{'type':'feet'},
            'shaft'  :{'type':'shaft'},
        },
    },
    'Venus': {
        'use_dual_drive': False,
        'pos': [0.0, 0.0, 1 * each_platter_z],
        'planetary_mult': 1,
        'scale_geneva': 0.85 * 0.625,
        'scale_planetary': 2.3 * 0.625,
        'rings_under_planets':True,
        'support_platter':'Luna',
        'rotor_azimuth':45,
        'scale_base_ring_radius':1.0,
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':29},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':14, 'inner_rail':0},
            'Pr'     :{'type':'spur',   'inout':'in', 'teeth':38, 'outer_rail':1},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':31},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':38, 'outer_rail':1},
            'Dr1'     :{'type':'spur',   'inout':'out', 'teeth':int(38 * 0.33)},
            'Dr2'     :{'type':'spur',   'inout':'out', 'teeth':int(38 * 0.33)},
            'Grotor' :{'type':'rotor'},
            'feet'   :{'type':'feet'},
            'shaft'  :{'type':'shaft'},
        },
    },
    'Mercury': {
        'use_dual_drive': False,
        'pos': [0.0, 0.0, 2 * each_platter_z],
        'planetary_mult': 1,
        'scale_geneva': 0.85 * 0.7,
        'scale_planetary': 2.0 * 0.7,
        'rings_under_planets':True,
        'support_platter':'Venus',
        'rotor_azimuth':60,
        'scale_base_ring_radius':1.0,
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':13},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':21, 'inner_rail':0},
            'Pr'     :{'type':'spur',   'inout':'in', 'teeth':39, 'outer_rail':1},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':29},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':23, 'outer_rail':1},
            'Dr1'     :{'type':'spur',   'inout':'out', 'teeth':int(23 * 0.33)},
            'Dr2'     :{'type':'spur',   'inout':'out', 'teeth':int(23 * 0.33)},
            'Grotor' :{'type':'rotor'},
            'feet'   :{'type':'feet'},
            'shaft'  :{'type':'shaft'},
        },
    },
}

def get_outer_tooth_shape(num_teeth, pitch_ref):
    original_shape = outer_tooth_shapes_p10[str(num_teeth)]
    shape = [0.1 * pitch_ref * np.array([-vert[0], vert[1], 0.0]) for vert in original_shape]
    return shape

def get_outer_tooth_inner_outer(num_teeth, pitch_ref):
    original_shape = outer_tooth_shapes_p10[str(num_teeth)]
    inner = pitch_ref * num_teeth / (2.0 * np.pi)
    inner_sq = inner * inner
    outer_sq = inner_sq
    scale = 0.1 * pitch_ref
    for vert in original_shape:
        r_sq = scale * vert[0] * scale * vert[0] + scale * vert[1] * scale * vert[1]
        if r_sq < inner_sq:
            inner_sq = r_sq
        if r_sq > outer_sq:
            outer_sq = r_sq
    inner = np.sqrt(inner_sq)
    outer = np.sqrt(outer_sq)
    return (inner, outer)


# def make_tooth_table(gear):
#     # https://en.wikipedia.org/wiki/Involute_gear
#     # Symbols from https://khkgears.net/new/gear_knowledge/gear_technical_reference/involute_gear_profile.html
#     num_teeth = gear['teeth']
#     pitch_ref = spur_tooth_pitch

#     p = pitch_ref
#     m = p / np.pi # Module is the unit size indicated in millimeter (mm)
#     ha = 1.0 * m # The distance between reference line and tooth tip
#     hf = 1.25 * m # The distance between reference line and tooth root
#     h = ha + hf # The distance between tooth tip and tooth root
#     hw = 2.0 * m # Depth of tooth meshed with the mating gear
#     c = 0.25 * m # The distance (clearance) between tooth root and the tooth tip of mating gear
#     rhof = 0.38 * m # The radius of curvature between tooth surface and the tooth root
#     rc = p * num_teeth / (2.0 * np.pi) # circular (reference) radius
#     ri = rc - hf # inner (root) radius
#     ro = rc + ha # outer (tip) radius

#     gear['specs'] = {
#                      'pitch_ref':pitch_ref,
#                      'radius_ref':rc,
#                      'radius_inner':ri,
#                      'radius_outer':ro,
#                      }
#     face_segments = 10
#     pitch_theta = 2.0 * np.pi / num_teeth
#     print('num_teeth', num_teeth)
#     print('radius_ref', rc, 'radius_inner', ri, 'radius_outer', ro)
#     print('pitch_theta', pitch_theta)
#     print('')
#     in_theta = 0.0
#     ref_theta = get_involute_theta(rc, ri)
#     out_theta = get_involute_theta(ro, ri)
#     total_tooth_theta = 2.0 * np.pi / num_teeth
#     inner_flat_theta = 0.5 * total_tooth_theta - 2.0 * ref_theta
#     outer_flat_theta = 0.5 * total_tooth_theta - 2.0 * (out_theta - ref_theta)
#     print('ref_theta',ref_theta)
#     print('out_theta',out_theta)
#     print('inner_flat_theta',inner_flat_theta)
#     print('outer_flat_theta',outer_flat_theta)
#     print('total_tooth_theta',total_tooth_theta)

#     tooth_curve_segments = 10
#     # tc_radius1 = []
#     # tc_theta1  = []
#     # for i in range(tooth_curve_segments):
#     #     t = float(i) / tooth_curve_segments
#     #     r = (rc - ri) * t + ri
#     #     theta = get_involute_theta(r, ri)
#     #     tc_radius1.append(r)
#     #     tc_theta1.append(theta)
#     # tc_radius2 = []
#     # tc_theta2  = []
#     # for i in range(tooth_curve_segments):
#     #     t = float(i) / tooth_curve_segments
#     #     r = (ro - rc) * t + rc
#     #     theta = get_involute_theta(r, ri)
#     #     tc_radius2.append(r)
#     #     tc_theta2.append(theta)
#     result_theta = []
#     result_radius = []
#     # Inner flat
#     theta_start = 0.0
#     for i in range(tooth_curve_segments):
#         t = float(i) / tooth_curve_segments
#         theta = inner_flat_theta * t
#         result_theta.append(theta + theta_start)
#         result_radius.append(ri)
#     theta_start += inner_flat_theta
#     # Tooth rise
#     for i in range(tooth_curve_segments):
#         t = float(i) / tooth_curve_segments
#         r = (ro - ri) * t + ri
#         theta = get_involute_theta(r, ri)
#         result_theta.append(theta + theta_start)
#         result_radius.append(r)
#     theta_start += out_theta
#     # Outer flat
#     for i in range(tooth_curve_segments):
#         t = float(i) / tooth_curve_segments
#         theta = outer_flat_theta * t
#         result_theta.append(theta + theta_start)
#         result_radius.append(ro)
#     theta_start += outer_flat_theta
#     # Tooth fall
#     for i in range(tooth_curve_segments):
#         t = float(i) / tooth_curve_segments
#         r = (ro - ri) * (1.0 - t) + ri
#         theta = get_involute_theta(r, ri)
#         result_theta.append((out_theta - theta) + theta_start)
#         result_radius.append(r)


#     return[result_theta, result_radius]

# def get_involute_theta(radius, root_radius):
#     alpha = np.arccos(root_radius / radius)
#     involute_alpha = np.tan(alpha) - alpha
#     return np.arcsin(involute_alpha)

####### This is *almost* right, but
####### the correect version can be found here: http://hessmer.org/gears/InvoluteSpurGearBuilder.html
def build_one_spur(gear):
    num_teeth = gear['teeth']
    specs = gear['specs']
    tooth_table2 = get_outer_tooth_shape(num_teeth, specs['pitch_ref'])
    num_verts = 2 * num_teeth * len(tooth_table2)
    verts = np.ndarray((num_verts, 3), dtype=np.float64)

    if gear.get('inout') == 'in':
        outer_mult = -1.0
        inner_vert_offset = 0
        outer_vert_offset = 1
    else:
        outer_mult = 1.0
        inner_vert_offset = 1
        outer_vert_offset = 0

    # Outer verts
    for tooth in range(num_teeth):
        tooth_t = float(tooth) / float(num_teeth)
        tooth_start_theta = 2.0 * np.pi * tooth_t
        tooth_vert_index = 2 * tooth * len(tooth_table2)
        sval = np.sin(tooth_start_theta)
        cval = np.cos(tooth_start_theta)
        for i in range(len(tooth_table2)):
            tx = tooth_table2[i][0]
            ty = tooth_table2[i][1]
            outer_vert_index = tooth_vert_index + i*2 + outer_vert_offset
            verts[outer_vert_index][0] = tx * cval + ty * sval
            verts[outer_vert_index][1] = ty * cval - tx * sval
            verts[outer_vert_index][2] = 0.0

    # Inner verts
    for tooth in range(num_teeth):
        tooth_vert_index = 2 * tooth * len(tooth_table2)
        for i in range(len(tooth_table2)):
            outer_vert_index = tooth_vert_index + i*2 + outer_vert_offset
            inner_vert_index = tooth_vert_index + i*2 + inner_vert_offset
            prev_outer_vert = verts[(outer_vert_index - 2 + num_verts) % num_verts]
            this_outer_vert = verts[outer_vert_index]
            next_outer_vert = verts[(outer_vert_index + 2) % num_verts]
            d1 = normalized(prev_outer_vert - this_outer_vert)
            d2 = normalized(next_outer_vert - this_outer_vert)
            c = np.cross(d1, d2)
            inner_vert_dir = outer_mult * np.sign(c[2]) * normalized(d1 + d2)
            verts[inner_vert_index] = this_outer_vert + tooth_rim_thickness * inner_vert_dir
#            verts[inner_vert_index] = this_outer_vert * 0.95
    gear['verts'] = [verts]
    gear['verts_z'] = [[-0.5 * spur_teeth_thickness, 0.5 * spur_teeth_thickness]]
    outer_rail_width = gear.get('outer_rail', None)
    if outer_rail_width is not None:
        rail_in = gear['specs']['radius_outer'] + 0.5 * thinnest_material_wall
        rail_out1 = rail_in + outer_rail_width * thinnest_material_wall
        rail_out2 = rail_in + 2 * thinnest_material_wall
        gear['verts'] += [make_cylinder_verts(rail_in, rail_out1)]
        gear['verts_z'] += [[-0.5 * spur_teeth_thickness - slide_buffer_dist, 0.5 * spur_teeth_thickness + slide_buffer_dist]]
        # gear['verts'] += [make_cylinder_verts(rail_in, rail_out2)]
        # gear['verts_z'] += [[-0.5 * spur_teeth_thickness, 0.5 * spur_teeth_thickness]]

    if gear['inout'] == 'out':
        if gear.get('inner_rail', 1):
            rail_out = gear['specs']['radius_inner'] - 0.5 * thinnest_material_wall
            rail_in = rail_out - 2 * thinnest_material_wall
            gear['axle_radius'] = rail_in
            gear['verts'] += [make_cylinder_verts(rail_in, rail_out)]
            gear['verts_z'] += [[-0.5 * spur_teeth_thickness - slide_buffer_dist, 0.5 * spur_teeth_thickness + slide_buffer_dist]]

    if gear.get('ring_under', None) is not None:
        rail_out = gear['specs']['radius_outer'] + 0.5 * thinnest_material_wall
        rail_in = rail_out - thinnest_material_wall
        rail_top = -0.5 * spur_teeth_thickness + 0.5 * thinnest_material_wall
        rail_bottom = rail_top - thinnest_material_wall
        gear['verts'] += [make_cylinder_verts(rail_in, rail_out)]
        gear['verts_z'] += [[rail_bottom, rail_top]]

    return verts

def make_geneva_tooth_table(gear, rotor):
    num_teeth = gear['teeth']
    pitch_ref = gear['specs']['pitch_ref']
    r = gear['specs']['radius_ref']
    total_tooth_theta = 2.0 * np.pi / num_teeth
    num_segments = 20
    groove_bottom_radius = r - rotor['arm_length'] + rotor['outset']
#    groove_bottom_radius = r - 0.1
    groove_half_width = rotor['pin_radius']
    groove_half_theta_inner = np.arctan(groove_half_width/groove_bottom_radius)
    groove_half_theta_outer = np.arctan(groove_half_width/r)
    result_radius = []
    result_theta = []

    # Curved groove bottom
    small_curve_segs = 12
    for seg in range(1, small_curve_segs):
        t = float(seg) / float(small_curve_segs)
        theta = np.radians(t * 90)
        x = rotor['pin_radius'] * np.sin(theta)
        y = groove_bottom_radius - rotor['pin_radius'] * np.cos(theta)
        result_radius += [y]
        result_theta += [np.arctan(x/y)]

    # Groove 1
    result_radius += [groove_bottom_radius]
    result_theta += [groove_half_theta_inner]
    result_radius += [r]
    result_theta += [groove_half_theta_outer]

    # rotor hub span
    span_segs = 30
    d = rotor['outset'] + r
    for sseg in range(1, span_segs):
        sst = float(sseg) / float(span_segs)
        seg_theta = sst * (0.5 * total_tooth_theta - groove_half_theta_outer) + groove_half_theta_outer
        seg_radius = r
        h = d * np.abs(np.sin(0.5 * total_tooth_theta - seg_theta))
        if h < rotor['hub_radius']:
            new_r = d - np.sqrt(rotor['hub_radius'] * rotor['hub_radius'] - h * h)
            if new_r < seg_radius:
                seg_radius = new_r
        result_radius += [seg_radius]
        result_theta += [seg_theta]

    # Replay everything in reverse
    result_theta += [total_tooth_theta - t for t in result_theta[::-1]]
    result_radius += [r for r in result_radius[::-1]]

    return [result_theta, result_radius]

def build_one_geneva(gear, rotor):
    tooth_table = make_geneva_tooth_table(gear, rotor)
    num_teeth = gear['teeth']
    num_verts = 2 * num_teeth * len(tooth_table[0])
    verts = np.ndarray((num_verts, 3), dtype=np.float64)
    # print(verts)
    # print(len(verts))
    # print(verts[0])
    # exit()
    specs = gear['specs']
    # Outer verts
    for tooth in range(num_teeth):
        tooth_t = float(tooth) / float(num_teeth)
        tooth_start_theta = 2.0 * np.pi * tooth_t
        tooth_vert_index = 2 * tooth * len(tooth_table[0])
        for i in range(len(tooth_table[0])):
            theta = tooth_start_theta + tooth_table[0][i]
            radius = tooth_table[1][i]
            # print(verts[tooth_vert_index + i*2])
            sval = np.sin(theta)
            cval = np.cos(theta)
            # print('theta', theta)
            # print('radius', radius)
            # print('sval', sval)
            outer_vert_index = tooth_vert_index + i*2
            verts[outer_vert_index][0] = radius * sval
            verts[outer_vert_index][1] = radius * cval
            verts[outer_vert_index][2] = 0.0
            # print(verts[tooth_vert_index + i*2])
            # print(verts[tooth_vert_index + i*2+1])
    # Inner verts
    for tooth in range(num_teeth):
        tooth_vert_index = 2 * tooth * len(tooth_table[0])
        for i in range(len(tooth_table[0])):
            outer_vert_index = tooth_vert_index + i*2
            inner_vert_index = outer_vert_index + 1
            prev_outer_vert = verts[(outer_vert_index - 2 + num_verts) % num_verts]
            this_outer_vert = verts[outer_vert_index]
            next_outer_vert = verts[(outer_vert_index + 2) % num_verts]
            d1 = normalized(prev_outer_vert - this_outer_vert)
            d2 = normalized(next_outer_vert - this_outer_vert)
            c = np.cross(d1, d2)
            inner_vert_dir = np.sign(c[2]) * normalized(d1 + d2)
            verts[inner_vert_index] = this_outer_vert + tooth_rim_thickness * inner_vert_dir
    gear['verts'] = [verts]
    gear['verts_z'] = [[-0.5 * thinnest_material_wall, 0.5 * thinnest_material_wall]]
    return verts

def configure_rotor(rotor, geneva):
    print('geneva[\'specs\'][\'pitch_ref\']',geneva['specs']['pitch_ref'])
    rotor['outset'] = geneva['specs']['pitch_ref'] * 0.2       # displacement from rim of G
    rotor['arm_length'] = np.sqrt((geneva['specs']['pitch_ref'] / 2.0) * (geneva['specs']['pitch_ref'] / 2.0) + rotor['outset'] * rotor['outset'])
    rotor['hub_radius'] = rotor['arm_length'] * 0.6   # radius of the hub
    rotor['pin_radius'] = rotor['arm_length'] * 0.2   # radius of the driver pin
    rotor['shaft_radius'] = rotor['hub_radius'] * 0.1   # radius of the shaft

def build_one_rotor(rotor, geneva):
    num_segments = 100
    num_verts = 2 * num_segments
    verts = np.ndarray((num_verts, 3), dtype=np.float64)
    pin_verts = np.ndarray((num_verts, 3), dtype=np.float64)
    disc_verts = np.ndarray((num_verts, 3), dtype=np.float64)
    # Outer verts
    for seg in range(num_segments):
        t = float(seg) / float(num_segments)
        theta = t * 2.0 * np.pi
        sval = np.sin(theta)
        cval = np.cos(theta)
        outer_radius = rotor['hub_radius'] - slide_buffer_dist
        inner_radius = outer_radius - thinnest_material_wall
        outer_vert_index = seg*2
        inner_vert_index = outer_vert_index + 1
        verts[outer_vert_index][0] = outer_radius * sval
        verts[outer_vert_index][1] = outer_radius * cval
        verts[outer_vert_index][2] = 0.0
        verts[inner_vert_index][0] = inner_radius * sval
        verts[inner_vert_index][1] = inner_radius * cval
        verts[inner_vert_index][2] = 0.0

        outer_pin_radius = rotor['pin_radius'] - slide_buffer_dist
        inner_pin_radius = rotor['pin_radius']  - thinnest_material_wall
        pin_verts[outer_vert_index][0] = rotor['arm_length'] + outer_pin_radius * sval
        pin_verts[outer_vert_index][1] = outer_pin_radius * cval
        pin_verts[outer_vert_index][2] = 0.0
        pin_verts[inner_vert_index][0] = rotor['arm_length'] + inner_pin_radius * sval
        pin_verts[inner_vert_index][1] = inner_pin_radius * cval
        pin_verts[inner_vert_index][2] = 0.0

        outer_disc_radius = 1.0 * rotor['pin_radius'] + rotor['arm_length']
        inner_disc_radius = outer_disc_radius - thinnest_material_wall
        disc_verts[outer_vert_index][0] = outer_disc_radius * sval
        disc_verts[outer_vert_index][1] = outer_disc_radius * cval
        disc_verts[outer_vert_index][2] = 0.0
        disc_verts[inner_vert_index][0] = inner_disc_radius * sval
        disc_verts[inner_vert_index][1] = inner_disc_radius * cval
        disc_verts[inner_vert_index][2] = 0.0

    rotor['verts'] = [verts, pin_verts, disc_verts]
    rotor['verts_z'] = [[-1.0 * thinnest_material_wall, 1.0 * spur_teeth_thickness],
                        [-1.0 * thinnest_material_wall, 1.0 * spur_teeth_thickness],
                        [-1.5 * thinnest_material_wall, -0.5 * thinnest_material_wall]]
    return verts

def normalized(v):
    return v / np.linalg.norm(v)
def length(v):
    return np.linalg.norm(v)

def write_one_quad(fp, v0, v1, v2, v3, n=None):
    if n is None:
        n = normalized(np.cross(v0 - v1, v2 - v1))
    fp.write('  facet normal {} {} {}\n'.format(n[0], n[1], n[2]))
    fp.write('    outer loop\n')
    fp.write('      vertex {} {} {}\n'.format(v0[0], v0[1], v0[2]))
    fp.write('      vertex {} {} {}\n'.format(v1[0], v1[1], v1[2]))
    fp.write('      vertex {} {} {}\n'.format(v2[0], v2[1], v2[2]))
    fp.write('    endloop\n')
    fp.write('  endfacet\n')
    fp.write('  facet normal {} {} {}\n'.format(n[0], n[1], n[2]))
    fp.write('    outer loop\n')
    fp.write('      vertex {} {} {}\n'.format(v2[0], v2[1], v2[2]))
    fp.write('      vertex {} {} {}\n'.format(v1[0], v1[1], v1[2]))
    fp.write('      vertex {} {} {}\n'.format(v3[0], v3[1], v3[2]))
    fp.write('    endloop\n')
    fp.write('  endfacet\n')

def write_one_tri(fp, v0, v1, v2, n=None):
    if n is None:
        n = normalized(np.cross(v0 - v1, v2 - v1))
    fp.write('  facet normal {} {} {}\n'.format(n[0], n[1], n[2]))
    fp.write('    outer loop\n')
    fp.write('      vertex {} {} {}\n'.format(v0[0], v0[1], v0[2]))
    fp.write('      vertex {} {} {}\n'.format(v1[0], v1[1], v1[2]))
    fp.write('      vertex {} {} {}\n'.format(v2[0], v2[1], v2[2]))
    fp.write('    endloop\n')
    fp.write('  endfacet\n')

def write_stl_ball_bearing(fp, collada_model, pos):
    global ball_bearing_faces
    if ball_bearing_faces is None:
        num_lat_segments = 6
        num_long_segments = 12 # must be even
        zflip = np.array([1.0, 1.0, -1.0])
        equator_verts = []
        rib_verts = []
        for seg in range(num_long_segments):
            t = seg / float(num_long_segments)
            theta = t * 2.0 * np.pi
            sval = np.sin(theta)
            cval = np.cos(theta)
            equator_verts.append(np.array([cval, sval, 1.0]))
        for seg in range(num_lat_segments):
            t = seg / float(num_lat_segments)
            theta = t * 0.5 * np.pi
            sval = np.sin(theta)
            cval = np.cos(theta)
            rib_verts.append(np.array([cval, cval, sval]))

        ball_bearing_faces = []
        for y in range(num_lat_segments - 1):
            rib0 = rib_verts[y + 0]
            rib1 = rib_verts[y + 1]
            for x in range(num_long_segments):
                eq0 = equator_verts[(x + 0) % num_long_segments]
                eq1 = equator_verts[(x + 1) % num_long_segments]
                q0 = eq0 * rib0 * ball_bearing_radius
                q1 = eq0 * rib1 * ball_bearing_radius
                q2 = eq1 * rib0 * ball_bearing_radius
                q3 = eq1 * rib1 * ball_bearing_radius
                ball_bearing_faces.append([q1, q0, q3, q2])
                ball_bearing_faces.append([q0*zflip, q1*zflip, q2*zflip, q3*zflip])
        # Now do the endcaps
        rib = rib_verts[num_lat_segments - 1]
        pole = np.array([0.0, 0.0, ball_bearing_radius])
        for x in range(num_long_segments):
            eq0 = equator_verts[(x + 0) % num_long_segments]
            eq1 = equator_verts[(x + 1) % num_long_segments]
            t0 = eq0 * rib * ball_bearing_radius
            t1 = eq1 * rib * ball_bearing_radius
            t2 = pole
            ball_bearing_faces.append([t0, t1, t2])
            ball_bearing_faces.append([t1*zflip, t0*zflip, t2*zflip])


    for face in ball_bearing_faces:
        if len(face) == 4:
            write_one_quad(fp, face[0]+pos, face[1]+pos, face[2]+pos, face[3]+pos)
        elif len(face) == 3:
            write_one_tri(fp, face[0]+pos, face[1]+pos, face[2]+pos)


def write_stl_tristrip_quads(fp, collada_model, verts, verts_z, pos=None, rot=None, closed=True):
    if pos is None:
        pos = np.array([0,0,0])
    if rot is None:
        rot = 0.0
    sval = np.sin(rot)
    cval = np.cos(rot)
    xverts = [np.array([v[0] * cval + v[1] * sval, v[1] * cval - v[0] * sval, v[2]]) for v in verts]
    xverts = xverts + pos
    num_verts = len(verts)
    num_quads = num_verts >> 1
    if not closed:
        num_quads -= 1
    z1 = np.array([0.0, 0.0, verts_z[1]])
    z2m1 = np.array([0.0, 0.0, verts_z[0] - verts_z[1]])
    n_up = np.array([0.0, 0.0, 1])
    n_down = np.array([0.0, 0.0, -1])
    for quad in range(num_quads):
        v0a = xverts[quad * 2] + z1
        v1a = xverts[quad * 2 + 1] + z1
        v2a = xverts[(quad * 2 + 2) % num_verts] + z1
        v3a = xverts[(quad * 2 + 3) % num_verts] + z1
        v0b = v0a + z2m1
        v1b = v1a + z2m1
        v2b = v2a + z2m1
        v3b = v3a + z2m1
        write_one_quad(fp, v0a, v1a, v2a, v3a, n=n_up)
        write_one_quad(fp, v1b, v0b, v3b, v2b, n=n_down)
        write_one_quad(fp, v0b, v0a, v2b, v2a)
        write_one_quad(fp, v1a, v1b, v3a, v3b)

def write_stl_platter(platter_name, collada_model, drive_gears_file_name, planet_gears_file_name, frame_file_name):
    platter = all_platters[platter_name]
    platter_pos = np.array(platter.get('pos', [0.0, 0.0, 0.0]))

    if collada_model is not None:
        material = collada_model.new_material((0,1,0),(0,1,0))

    # try:
    fp_drive_gears = open(drive_gears_file_name, 'w')
    fp_planet_gears = open(planet_gears_file_name, 'w')
    fp_frame = open(frame_file_name, 'w')

    for gear_name,gear in platter['gears'].items():
        print('  Writing gear {}:{}'.format(platter_name, gear_name))
        if gear_name in ['feet']:
            fp = fp_frame
        elif gear_name in ['G', 'Grotor', 'Da', 'Db']:
            fp = fp_drive_gears
        else:
            fp = fp_planet_gears
        fp.write('solid OpenSCAD_Model\n')
        verts = gear.get('verts', None)
        if verts is not None:
            for strip_index,strip_verts in enumerate(verts):
                pos = gear.get('pos', [0.0, 0.0, 0.0]) + platter_pos
                rot = gear.get('rot', 0.0)
                verts_z = gear['verts_z'][strip_index]
                write_stl_tristrip_quads(fp, collada_model, strip_verts, verts_z=verts_z, pos=pos, rot=rot, closed=True)
                if collada_model is not None:
                    collada_model.add_extruded_tristrip_quad_mesh(material, strip_verts, verts_z=verts_z, pos=pos, rot=rot, closed=True)
        bearings = gear.get('ball_bearings', None)
        if bearings is not None:
            for bearing_pos in bearings:
                write_stl_ball_bearing(fp, collada_model, pos=bearing_pos)

        fp.write('endsolid OpenSCAD_Model\n')
    fp_drive_gears.close()
    fp_planet_gears.close()
    fp_frame.close()
    # except:
    #     print('Unable to write file', file_name)

def make_platter_specs(platter_name):
    platter = all_platters[platter_name]
    for gear_name,gear in platter['gears'].items():
        if gear_name in ['Ps', 'Pr', 'Pp0', 'Pp1', 'Pp2', 'Pp3']:
            planetary_mult = platter.get('planetary_mult', 1)
            gear['teeth'] *= planetary_mult
        gear['specs'] = {}
        gear['rings'] = []
        if gear['type'] == 'geneva':
            gear['specs']['radius_ref'] = geneva_outer_radius
            gear['specs']['pitch_ref'] = 2.0 * np.pi * gear['specs']['radius_ref'] / gear['teeth']
        elif gear['type'] == 'spur':
            pitch_ref = spur_tooth_pitch
            num_teeth = gear['teeth']
            gear['specs']['pitch_ref'] = pitch_ref
            gear['specs']['radius_ref'] = pitch_ref * num_teeth / (2.0 * np.pi)
            inner,outer = get_outer_tooth_inner_outer(num_teeth, pitch_ref)
            gear['specs']['radius_inner'] = inner
            gear['specs']['radius_outer'] = outer

    # Find the suitable planet angles
    if 'Ps' in platter['gears'] and 'Pr' in platter['gears']:
        sun = platter['gears']['Ps']
        ring = platter['gears']['Pr']
        platter['planetary_nooks'] = []
        sun_steps = 2 * sun['teeth']
        ring_steps = 2 * ring['teeth']
        for sun_step in range(sun_steps):
            t = float(sun_step) / float(sun_steps)
            ring_step = t * ring_steps
            err = np.abs(ring_step - np.round(ring_step, 1))
            if err < 0.0001:
                platter['planetary_nooks'].append((t * 2.0 * np.pi, sun_step & 1))
#            print('   step {},{} err {}'.format(sun_step, ring_step, err))
        print('    planetary_nooks:',len(platter['planetary_nooks']))
        skip = 1
        num_planet_gears = len(platter['planetary_nooks'])
        skip = 1
        while num_planet_gears > 4:
            num_planet_gears >>= 1
            skip <<= 1
        for i in range(num_planet_gears):
            nook = platter['planetary_nooks'][i * skip]
            p = {'type':'spur', 'inout':'out', 'specs':{}}
            p['teeth'] = (ring['teeth'] - sun['teeth']) >> 1
            p['specs']['pitch_ref'] = sun['specs']['pitch_ref']
            p['specs']['radius_ref'] = p['specs']['pitch_ref'] * p['teeth'] / (2.0 * np.pi)
            p['p_placement'] = nook
            inner,outer = get_outer_tooth_inner_outer(p['teeth'], p['specs']['pitch_ref'])
            p['specs']['radius_inner'] = inner
            p['specs']['radius_outer'] = outer
            if platter.get('rings_under_planets', 0):
                p['ring_under'] = 1
            platter['gears']['Pp{}'.format(i)] = p


def make_cylinder_verts(inner_radius, outer_radius, center=None, num_segments=None):
    if num_segments is None:
        num_segments = 100
    num_verts = 2 * num_segments
    verts = np.ndarray((num_verts, 3), dtype=np.float64)

    for seg in range(num_segments):
        t = float(seg) / float(num_segments)
        theta = t * 2.0 * np.pi
        sval = np.sin(theta)
        cval = np.cos(theta)
        outer_vert_index = seg*2
        inner_vert_index = outer_vert_index + 1
        verts[outer_vert_index][0] = outer_radius * sval
        verts[outer_vert_index][1] = outer_radius * cval
        verts[outer_vert_index][2] = 0.0
        verts[inner_vert_index][0] = inner_radius * sval
        verts[inner_vert_index][1] = inner_radius * cval
        verts[inner_vert_index][2] = 0.0

        if center is not None:
            verts[outer_vert_index] += center
            verts[inner_vert_index] += center
    return verts

def add_ball_bearings(gear, count, center, placement_radius):
    bearings = gear.get('ball_bearings', [])
    for i in range(count):
        t = i / float(count)
        theta = 2.0 * np.pi * t;
        sval = np.sin(theta)
        cval = np.cos(theta)
        v = np.array([center[0] + cval * placement_radius,
                      center[1] + sval * placement_radius,
                      center[2]])
        bearings.append(v)
    gear['ball_bearings'] = bearings


def make_gear_hub_shaft(platter_name, gear, strut_bottom, strut_top):
    platter = all_platters[platter_name]

    strut_verts = make_cylinder_verts(strut_inner_radius, strut_outer_radius)
    strut_pos = np.array([gear['pos'][0], gear['pos'][1], 0.0])
    feet = platter['gears']['feet']
    # feet['verts'].append(strut_verts + strut_pos)
    # feet['verts_z'].append([strut_bottom, strut_top])

    ring_radius = platter['base_ring_radius']
    if platter_name == 'Luna':
        ring_radius = all_platters['Mars']['base_ring_radius']

    origin = np.array([0.0, 0.0, 0.0])


    ## make the axle
    axle1_radius = gear['axle_radius'] - slide_buffer_dist
    axle1_in = axle1_radius - thinnest_material_wall
    axle1_out = axle1_radius
    axle1_verts = make_cylinder_verts(axle1_in, axle1_out)
    axle1_pos = np.array([gear['pos'][0], gear['pos'][1], 0.0])
    axle1_top = gear['pos'][2] + 0.5 * tooth_rim_thickness
    axle1_bottom = axle1_top - tooth_rim_thickness - 0.5 * thinnest_material_wall
    shaft = platter['gears']['shaft']
    feet['verts'].append(axle1_verts + axle1_pos)
    feet['verts_z'].append([axle1_bottom, axle1_top])

    bearing_center = axle1_pos + (0.0, 0.0, 0.5 * (axle1_top + axle1_bottom) + platter['pos'][2])
    placement_radius = axle1_out - 0.5 * ball_bearing_radius
    bearing_count = int(2.0 * np.pi * placement_radius / (4 * ball_bearing_radius))
    add_ball_bearings(feet, count = bearing_count, center=bearing_center, placement_radius=placement_radius)

    # make the axle bottom-rest
    axle2_in = axle1_in
    axle2_out = 0.5 * gear['specs']['radius_ref'] + 0.5 * gear['specs']['radius_inner']
    if axle2_out < axle2_in + 2 * strut_outer_radius:
        axle2_out = axle2_in + 2 * strut_outer_radius
    axle2_top = axle1_top - tooth_rim_thickness - slide_buffer_dist
    gear['axle_rest_z'] = axle2_top
    axle2_bottom = axle2_top - thinnest_material_wall
    axle2_verts = make_cylinder_verts(axle2_in, axle2_out)
    feet['verts'].append(axle2_verts + axle1_pos)
    feet['verts_z'].append([axle2_bottom, axle2_top])

    bearing_center = axle1_pos + (0.0, 0.0, axle2_top - 0.5 * ball_bearing_radius + platter['pos'][2])
    placement_radius = gear['specs']['radius_inner'] - 2.0 * ball_bearing_radius
    bearing_count = int(2.0 * np.pi * placement_radius / (4 * ball_bearing_radius))
    add_ball_bearings(feet, count = bearing_count, center=bearing_center, placement_radius=placement_radius)

    # Struts
    if 0:
        strut_top = axle2_bottom + slide_buffer_dist
        mini_radius = axle2_in + 1 * strut_outer_radius
        ring_pts = find_circle_intersection_points(origin, ring_radius, gear['pos'], mini_radius)
        if len(ring_pts):
            feet['verts'].append(strut_verts + ring_pts[0])
            feet['verts_z'].append([strut_bottom, strut_top])
            if len(ring_pts) > 1:
                feet['verts'].append(strut_verts + ring_pts[1])
                feet['verts_z'].append([strut_bottom, strut_top])

def make_gear_hub_riser(platter_name, gear, strut_bottom, strut_top):
    platter = all_platters[platter_name]
    strut_verts = make_cylinder_verts(strut_inner_radius, strut_outer_radius)
    strut_pos = np.array([gear['pos'][0], gear['pos'][1], 0.0])
    feet = platter['gears']['feet']
    feet['verts'].append(strut_verts + strut_pos)
    feet['verts_z'].append([strut_bottom, strut_top])

def make_platter_supports(platter_name):
    platter = all_platters[platter_name]
    feet = platter['gears']['feet']
    shaft = platter['gears']['shaft']
    gearDa = platter['gears'].get('Da', None)
    gearPs = platter['gears'].get('Ps', None)
    foot_verts = make_cylinder_verts(strut_inner_radius, strut_outer_radius)
    r2 = 1.0 / np.sqrt(2.0)
    if gearDa:
        shaft_gear = gearDa
    else:
        shaft_gear = gearPs
    fd = r2 * (shaft_gear['specs']['radius_outer'] + 1 * strut_outer_radius + thinnest_material_wall)
    sd = r2 * (shaft_gear['specs']['radius_ref'] - 1.5 * strut_outer_radius)
    f2d = r2 * (gearPs['specs']['radius_inner'] - 1 * strut_outer_radius - thinnest_material_wall - slide_buffer_dist)
    foot_offsets = [np.array([-1, -1, 0.0]),
                    np.array([ 1, -1, 0.0]),
                    np.array([-1,  1, 0.0]),
                    np.array([ 1,  1, 0.0])]
    feet['verts'] = []
    feet['verts_z'] = []
    shaft['verts'] = []
    shaft['verts_z'] = []

    support_platter_name = platter['support_platter']
    support_platter = all_platters.get(support_platter_name, None)
    base_clearance = 2.0 * spur_teeth_thickness
    base_z = platter['pos'][2] - base_clearance

    platter['base_ring_radius'] = length(platter['gears']['Pp0']['pos']) * platter.get('scale_base_ring_radius', 1.0)

    for gear_name,gear in platter['gears'].items():
        strut_bottom = platter['base_z']
        if gear_name in ['Db']:
            strut_top = gear['pos'][2]
            make_gear_hub_shaft(platter_name, gear, strut_bottom, strut_top)
        if gear_name in ['Pp0', 'Pp1', 'Pp2', 'Pp3']:
            strut_top = gear['pos'][2] + 8
            make_gear_hub_shaft(platter_name, gear, strut_bottom, strut_top)
        if platter_name in ['Drive'] and gear_name in ['Ps']:
            strut_top = gear['pos'][2] + 8
            make_gear_hub_shaft(platter_name, gear, strut_bottom, strut_top)

    if 0:
        for gear_name,gear in platter['gears'].items():
            if gear['type'] == 'feet':
                for foot_index in range(4):
                    feet['verts'].append(foot_verts + fd * foot_offsets[foot_index])
                    feet['verts_z'].append([shaft_gear['pos'][2] - 1.0, shaft_gear['pos'][2] + 1.0])
                    feet['verts'].append(foot_verts + f2d * foot_offsets[foot_index])
                    feet['verts_z'].append([gearPs['pos'][2] - 1.0, gearPs['pos'][2] + 1.0])
                    shaft['verts'].append(foot_verts + sd * foot_offsets[foot_index])
                    shaft['verts_z'].append([shaft_gear['pos'][2] - 1.0, shaft_gear['pos'][2] + 1.0])

def set_platter_base(platter_name):
    platter = all_platters[platter_name]
    if platter_name == 'Drive':
        platter['base_z'] = -1.5 * thinnest_material_wall
    elif platter_name == 'Mars':
        platter['base_z'] = all_platters['Drive']['base_z'] - each_platter_z
    elif platter_name == 'Luna':
        platter['base_z'] = all_platters['Drive']['base_z'] - 2 * each_platter_z
    else:
        platter['base_z'] = -2.0


def do_platter_adjustments(platter_name):
    platter = all_platters[platter_name]
    if platter_name == 'Drive':
        platter['gears']['G']['pos'][2] = platter['gears']['Ps']['pos'][2]

    if platter_name in ['Drive']:
        riser_height = None
        riser_plan = None
        riser_rings = None
        gears = platter['gears']
        rad_in = gears['Ps']['axle_radius']
        rad_out = platter['base_ring_radius'] + gears['Pp0']['axle_radius']
        geneva_ring_radius = 0.25 * gears['Pr']['outer_rail'] + gears['Pr']['specs']['radius_outer']
        center_rings = [[rad_out, thinnest_material_wall], [geneva_ring_radius, 1.5*thinnest_material_wall]]
        connect_out = all_platters['Mars']['base_ring_radius'] + all_platters['Mars']['gears']['Pp0']['axle_radius']
        connect_rings = [[[rad_in, connect_out],[0, 90, 180, 270]]]
        ring_radius = None
        ring_base_z = platter['base_z']
        make_base_ring(platter_name, ring_radius, ring_base_z,
                       center_rings=center_rings, connect_rings=connect_rings,
                       riser_height=riser_height, riser_rings=riser_rings)
    if platter_name in ['Mars']:
        rad_in = all_platters['Drive']['base_ring_radius'] + all_platters['Drive']['gears']['Pp0']['axle_radius']
        rad_out = platter['base_ring_radius'] + platter['gears']['Pp0']['axle_radius']
        rad_mid = platter['base_ring_radius'] - platter['gears']['Pp0']['axle_radius']
        rad_out2 = 0.5 * rad_out + 0.5 * rad_mid
        # ring_radius = platter['base_ring_radius']
        ring_radius = None
        ring_base_z = platter['base_z']
        riser_height = 1.0
        riser_height_Pp = platter['gears']['Pp0']['axle_rest_z'] - 1.0 * thinnest_material_wall - ring_base_z
        riser_height_Db = platter['gears']['Db']['axle_rest_z'] - 1.0 * thinnest_material_wall - ring_base_z
        riser_plan = [
                     [0.0-5.0,  riser_height_Pp], 
                     [0.0+5.0,  riser_height_Pp],
                     [0.0+20.0, 0.0],
                     [180.0-20.0, 0.0], 
                     [180.0-5.0,  riser_height_Pp], 
                     [180.0+5.0,  riser_height_Pp],
                     [180.0+20.0, 0.0],
                     [360.0-20.0, 0.0], 
                     [360.0-5.0,  riser_height_Pp], 
                     [360.0,      riser_height_Pp],
                     ]
        riser_planA = [
                     [0.0-5.0,  riser_height_Pp], 
                     [0.0+5.0, 0.0],
                     [180.0-15.0, 0.0], 
                     [180.0-5.0,  riser_height_Pp], 
                     [180.0+5.0, 0.0], 
                     [360.0-15.0, 0.0], 
                     [360.0-5.0,  riser_height_Pp], 
                     [360.0+5.0,   0.0],
                     ]
        riser_planB = [
                     [0.0-5.0,  0.0], 
                     [0.0+5.0,  riser_height_Pp], 
                     [0.0+15.0, 0.0],
                     [180.0-5.0, 0.0], 
                     [180.0+5.0,  riser_height_Pp], 
                     [180.0+15.0, 0.0], 
                     [360.0-5.0, 0.0], 
                     [360.0+5.0,  riser_height_Pp], 
                     ]
        riser_planAA = [
                     [0.0-10.0,  riser_height_Pp], 
                     [0.0+0.0, 0.0],
                     [180.0-20.0, 0.0], 
                     [180.0-10.0,  riser_height_Pp], 
                     [180.0+0.0, 0.0], 
                     [360.0-20.0, 0.0], 
                     [360.0-10.0,  riser_height_Pp], 
                     [360.0+0.0,   0.0],
                     ]
        riser_planBB = [
                     [0.0-0.0,  0.0], 
                     [0.0+10.0,  riser_height_Pp], 
                     [0.0+20.0, 0.0],
                     [180.0-0.0, 0.0], 
                     [180.0+10.0,  riser_height_Pp], 
                     [180.0+20.0, 0.0], 
                     [360.0-0.0, 0.0], 
                     [360.0+10.0,  riser_height_Pp], 
                     ]
        riser_planC = [
                     [0.0,  0.0], 
                     [270.0-31.0+2.0-10.0,  0.0], 
                     [270.0-31.0+2.0,  riser_height_Db], 
                     [270.0-31.0+2.0+10.0,  0.0], 
                     [270.0+31.0+2.0-10.0,  0.0], 
                     [270.0+31.0+2.0,  riser_height_Db], 
                     [270.0+31.0+2.0+10.0,  0.0], 
                     [360.0,  0.0], 
                     ]
        riser_planD = [
                     [0.0,  0.0], 
                     [270.0-31.0-2.0-10.0,  0.0], 
                     [270.0-31.0-2.0,  riser_height_Db], 
                     [270.0-31.0-2.0+10.0,  0.0], 
                     [270.0+31.0-2.0-10.0,  0.0], 
                     [270.0+31.0-2.0,  riser_height_Db], 
                     [270.0+31.0-2.0+10.0,  0.0], 
                     [360.0,  0.0], 
                     ]
        riser_planE = [
                     [0.0,  0.0], 
                     [270.0-37.0+2.0-10.0,  0.0], 
                     [270.0-37.0+2.0,  riser_height_Db], 
                     [270.0-37.0+2.0+10.0,  0.0], 
                     [270.0+37.0+2.0-10.0,  0.0], 
                     [270.0+37.0+2.0,  riser_height_Db], 
                     [270.0+37.0+2.0+10.0,  0.0], 
                     [360.0,  0.0], 
                     ]
        riser_planF = [
                     [0.0,  0.0], 
                     [270.0-37.0-2.0-10.0,  0.0], 
                     [270.0-37.0-2.0,  riser_height_Db], 
                     [270.0-37.0-2.0+10.0,  0.0], 
                     [270.0+37.0-2.0-10.0,  0.0], 
                     [270.0+37.0-2.0,  riser_height_Db], 
                     [270.0+37.0-2.0+10.0,  0.0], 
                     [360.0,  0.0], 
                     ]
        center_rings = [[rad_out, thinnest_material_wall]]
#        connect_rings = [[[rad_in, rad_out],[0, 90, 180, 270]]]
        connect_rings = None
        riser_rings = [
                       [rad_out, riser_planA],
                       [rad_out, riser_planB],
                       [rad_out2, riser_planAA],
                       [rad_out2, riser_planBB],
                       [rad_out, riser_planC],
                       [rad_out, riser_planD],
                       [rad_out2, riser_planE],
                       [rad_out2, riser_planF],
                       ]
        make_base_ring(platter_name, ring_radius, ring_base_z,
                       center_rings=center_rings, connect_rings=connect_rings,
                       riser_height=riser_height, riser_rings=riser_rings)
    elif platter_name in ['Venus', 'Mercury']:
        riser_height = None
        riser_plan = None
        center_rings = None
        connect_rings = None
        riser_rings = None
        ring_radius = platter['base_ring_radius']
        ring_base_z = platter['base_z']
        make_base_ring(platter_name, ring_radius, ring_base_z,
                       center_rings=center_rings, connect_rings=connect_rings,
                       riser_height=riser_height, riser_rings=riser_rings)

def find_circle_intersection_points(c1, r1, c2, r2):
    d = length(c2 - c1)
    if d >= r1 + r2:
        return []
    a = (r1*r1 - r2*r2 + d*d) / (2*d)
    cross_pt = c1 + (c2 - c1) * a / d
    if r1*r1 < a*a:
        return []
    h = np.sqrt(r1*r1 - a*a)
    dx =  h * (c2[1] - c1[1]) / d
    dy = -h * (c2[0] - c1[0]) / d
    dxy = np.array([dx, dy, 0.0])
    return [cross_pt + dxy, cross_pt - dxy]

def make_base_ring(platter_name, ring_radius, ring_base_z,
                   center_rings=None, connect_rings=None, 
                   riser_height=None, riser_rings=None):
    platter = all_platters[platter_name]
    ring_top = ring_base_z + 1.0 * thinnest_material_wall
    ring_bottom = ring_top - thinnest_material_wall
    feet = platter['gears']['feet']

    span_dist = strut_outer_radius * 6

    if center_rings is not None:
        for ring in center_rings:
            this_ring_top = ring_bottom + ring[1]
            ring_out = ring[0] + 0.5 * thinnest_material_wall
            ring_in = ring_out - thinnest_material_wall
            ring_verts = make_cylinder_verts(ring_in, ring_out)
            feet['verts'].append(ring_verts)
            feet['verts_z'].append([ring_bottom, this_ring_top])

    if connect_rings is not None:
        for ring in connect_rings:
            connect_in = ring[0][0]
            connect_out = ring[0][1]
            thetas = ring[1]
            avg_rad = 0.5 * (connect_in + connect_out)
            ring_in = avg_rad - 0.5 * thinnest_material_wall
            ring_out = ring_in + thinnest_material_wall
            ring_offset_dist = avg_rad - connect_in
            ring_verts = make_cylinder_verts(ring_in, ring_out)
            for theta in thetas:
                sval = np.sin(np.radians(theta))
                cval = np.cos(np.radians(theta))
                ring_offset = np.array([ring_offset_dist * sval, ring_offset_dist * cval , 0.0])
                feet['verts'].append(ring_verts + ring_offset)
                feet['verts_z'].append([ring_bottom, ring_top])

    if riser_rings is not None:
        for ring in riser_rings:
            radius = ring[0]
            riser_plan = ring[1]
            ring_in = radius - 0.5 * thinnest_material_wall
            ring_out = radius + 0.5 * thinnest_material_wall
            ring_verts = make_cylinder_verts(ring_in, ring_out, num_segments=500)
            riser_verts = []
            for i,v in enumerate(ring_verts):
                th = 360.0 * float(i) / float(len(ring_verts))
                for plan_index,rp2 in enumerate(riser_plan):
                    if rp2[0] > th:
                        break
                rp1 = riser_plan[plan_index - 1]
                t = (th - rp1[0]) / (rp2[0] - rp1[0])
                if 0:
                    # use sin^2 to mellow it out
                    rise_t = np.sin(t * 0.5 * np.pi)
                    rise_t *= rise_t
                if 1:
                    # use sin to mellow it out
                    rise_t = np.sin((t - 0.5) * 1.0 * np.pi)
                    rise_t = 0.5 * (rise_t + 1.0)
                rise = rp1[1] + rise_t * (rp2[1] - rp1[1])

                riser_verts.append(v + [0.0, 0.0, riser_height * rise])

            feet['verts'].append(riser_verts)
            feet['verts_z'].append([ring_bottom, ring_top])

    ###################################
    ## This section is old and should be deleted
    if ring_radius is not None:
        # Outer ring
        ring_out1 = ring_radius + 0.5 * span_dist
        ring_in1 = ring_out1 - thinnest_material_wall
        ring_verts = make_cylinder_verts(ring_in1, ring_out1)
        feet['verts'].append(ring_verts)
        feet['verts_z'].append([ring_bottom, ring_top])
        # Outer riser
        # if riser_plan is not None:
        #     ring_verts = make_cylinder_verts(ring_in1, ring_out1, num_segments=500)
        #     riser_verts = []
        #     for i,v in enumerate(ring_verts):
        #         th = 360.0 * float(i) / float(len(ring_verts))
        #         for plan_index,rp2 in enumerate(riser_plan):
        #             if rp2[0] > th:
        #                 break
        #         rp1 = riser_plan[plan_index - 1]
        #         t = (th - rp1[0]) / (rp2[0] - rp1[0])
        #         rise = rp1[1] + t * (rp2[1] - rp1[1])
        #         if 0:
        #             # use sin^2 to mellow it out
        #             rise = np.sin(rise * 0.5 * np.pi)
        #             rise *= rise
        #         if 1:
        #             # use sin to mellow it out
        #             rise = np.sin((rise - 0.5) * 1.0 * np.pi)
        #             rise = 0.5 * (rise + 1.0)

        #         riser_verts.append(v + [0.0, 0.0, riser_height * rise])

        #     feet['verts'].append(riser_verts)
        #     feet['verts_z'].append([ring_bottom, ring_top])
        #return

        #inner ring
        ring_in2 = ring_radius - 0.5 * span_dist
        ring_out2 = ring_in2 + thinnest_material_wall
        ring_verts = make_cylinder_verts(ring_in2, ring_out2)
        feet['verts'].append(ring_verts)
        feet['verts_z'].append([ring_bottom, ring_top])
        #cross rings
        ring_in3 = ring_radius - 0.5 * thinnest_material_wall
        ring_out3 = ring_in3 + thinnest_material_wall
        ring_verts = make_cylinder_verts(ring_in3, ring_out3)
        ring_offset_dist = 0.5 * span_dist - 0.5 * thinnest_material_wall
        thetas = [0, 180, 90, 270]
        for theta in thetas:
            sval = np.sin(np.radians(theta))
            cval = np.cos(np.radians(theta))
            ring_offset = np.array([ring_offset_dist * sval, ring_offset_dist * cval , 0.0])
            feet['verts'].append(ring_verts + ring_offset)
            feet['verts_z'].append([ring_bottom, ring_top])


def tooth_theta(gear):
    return 2.0 * np.pi / gear['teeth']

def scale_gear_specs(gear, scale):
    if gear is None:
        return
    specs = gear['specs']
    if 'pitch_ref' in specs:
        specs['pitch_ref'] *= scale
    if 'radius_ref' in specs:
        specs['radius_ref'] *= scale
    if 'radius_inner' in specs:
        specs['radius_inner'] *= scale
    if 'radius_outer' in specs:
        specs['radius_outer'] *= scale

def build_one_platter(platter_name):
    platter = all_platters[platter_name]
    dual_drive = platter.get('use_dual_drive', False)
    gearG  = platter['gears'].get('G', None)
    gearPr = platter['gears'].get('Pr', None)
    gearPp = []
    for i in range(4):
        p = platter['gears'].get('Pp{}'.format(i), None)
        if p:
            gearPp.append(p)
    gearPs = platter['gears'].get('Ps', None)
    gearDa = platter['gears'].get('Da', None)
    gearDb = platter['gears'].get('Db', None)
    if not dual_drive:
        if 'Dr1' in platter['gears']:
            platter['gears'].pop('Dr1')
            platter['gears'].pop('Dr2')
    gearDr1 = platter['gears'].get('Dr1', None)
    gearDr2 = platter['gears'].get('Dr2', None)
    rotor = platter['gears']['Grotor']

    rotor_azimuth = platter.get('rotor_azimuth', 0.0)

    scale_planetary = platter.get('scale_planetary', None)
    if scale_planetary:
        scale_gear_specs(gearPr, scale_planetary)
        scale_gear_specs(gearPs, scale_planetary)
        for p in gearPp:
            scale_gear_specs(p, scale_planetary)

    scale_geneva = platter.get('scale_geneva', None)
    if scale_geneva:
        print('  AA ',scale_geneva)
        print('  A gearG[\'specs\'][\'pitch_ref\']',gearG['specs']['pitch_ref'])
        scale_gear_specs(gearG, scale_geneva)
        print('  A gearG[\'specs\'][\'pitch_ref\']',gearG['specs']['pitch_ref'])

    configure_rotor(rotor, gearG)

    if gearDa:
        if dual_drive:
            dadb_current = gearDa['specs']['radius_ref'] + gearDb['specs']['radius_ref'] + gearDr1['specs']['radius_ref'] + gearDr2['specs']['radius_ref']
        else:
            dadb_current = gearDa['specs']['radius_ref'] + gearDb['specs']['radius_ref']
        dadb_desired = gearG['specs']['radius_ref'] + rotor['outset']
        dadb_fix = dadb_desired / dadb_current
        scale_gear_specs(gearDa, dadb_fix)
        scale_gear_specs(gearDb, dadb_fix)
        if dual_drive:
            scale_gear_specs(gearDr1, dadb_fix)
            scale_gear_specs(gearDr2, dadb_fix)

        match_post_sizes = False
        if match_post_sizes:
            post1 = gearDa['specs']['radius_outer'] + 1 * strut_outer_radius + thinnest_material_wall
            post2 = gearPs['specs']['radius_inner'] - 0 * strut_outer_radius - 0 * thinnest_material_wall
            ring_fix = post1 / post2
            scale_gear_specs(gearPr, ring_fix)
            scale_gear_specs(gearPs, ring_fix)
            for p in gearPp:
                scale_gear_specs(p, ring_fix)

    match_ring_to_G = False
    if match_ring_to_G:
        if gearPr:
            ring_current = gearPr['specs']['radius_ref']
        ring_desired = planet_outer_ref_radius
        ring_fix = ring_desired / ring_current
        scale_gear_specs(gearPr, ring_fix)
        scale_gear_specs(gearPs, ring_fix)
        for p in gearPp:
            scale_gear_specs(p, ring_fix)

    for gear_name,gear in platter['gears'].items():
        gtype = gear['type']
        if gtype == 'geneva':
            build_one_geneva(gear, rotor)
        elif gtype == 'spur':
            build_one_spur(gear)
        elif gtype == 'rotor':
            build_one_rotor(gear, gearG)
    geneva_z = 0.0
    planetary_z = geneva_z + 0.5 * thinnest_material_wall + 0.5 * spur_teeth_thickness
    if dual_drive:
        mini_driver_z = geneva_z - 2.0 * thinnest_material_wall - 0.5 * spur_teeth_thickness
        mini_driver2_z = geneva_z - 1.5 * thinnest_material_wall - 0.5 * spur_teeth_thickness
        driver_z = mini_driver_z - 1.0 * spur_teeth_thickness - 0.0 * thinnest_material_wall
    else:
        driver_z = geneva_z - slide_buffer_dist - 0.5 * thinnest_material_wall - 0.5 * spur_teeth_thickness
    gearG['pos'] = np.array([0.0, 0.0, 0.0])
    gearG['rot'] = 0.5 * tooth_theta(gearG) + np.radians(rotor_azimuth)
    gearPs['pos'] = np.array([0.0, 0.0, planetary_z])
    if gearPr:
        gearPr['pos'] = np.array([0.0, 0.0, planetary_z])
        # gearPr['rot'] = 0.5 * tooth_theta(gearPr)
        gearPr['rot'] = 0.0
        if (gearPp[0]['teeth'] & 1) == 0:
            gearPr['rot'] = 0.5 * tooth_theta(gearPr)
    for pi,p in enumerate(gearPp):
        disp = gearPs['specs']['radius_ref'] + p['specs']['radius_ref']
        theta = p['p_placement'][0]
        # theta = 0
        rot = theta * float(gearPs['teeth']) / float(p['teeth'])
        if (p['teeth'] & 1) == 0:
            rot += 0.5 * tooth_theta(p)
        # hack special cases
        if platter_name == 'Mars' and pi == 1:
            rot += 0.5 * tooth_theta(p)
        # rot = 0
        sval = np.sin(theta)
        cval = np.cos(theta)
        p['pos'] = np.array([disp * sval, disp * cval, planetary_z])
        p['rot'] = rot

    if gearDa:
        bdist = gearDa['specs']['radius_ref'] + gearDb['specs']['radius_ref']
        gearDa['pos'] = np.array([0.0, 0.0, driver_z])
        gearDb['pos'] = np.array([0.0, bdist, driver_z])
        gearDb['rot'] = 0.0
        if (gearDb['teeth'] & 1) == 0:
            gearDb['rot'] += 0.5 * tooth_theta(gearDb)
        gearDb['rot'] += np.radians(rotor_azimuth) * 2.0 * float(gearDa['teeth']) / float(gearDb['teeth'])
        sval = np.sin(np.radians(rotor_azimuth))
        cval = np.cos(np.radians(rotor_azimuth))
        gearDb['pos'][0] = bdist * sval;
        gearDb['pos'][1] = bdist * cval;

        if dual_drive:
            gearDr1['pos'] = np.array([gearDb['pos'][0], gearDb['pos'][1], mini_driver_z])
            gearDr2['pos'] = np.array([0.0, gearDr1['pos'][1] + gearDr1['specs']['radius_ref'] + gearDr2['specs']['radius_ref'], mini_driver2_z])
            gearDr2['rot'] = 0.5 * tooth_theta(gearDr2)
    rotor_dist = gearG['specs']['radius_ref'] + rotor['outset']
    sval = np.sin(np.radians(rotor_azimuth))
    cval = np.cos(np.radians(rotor_azimuth))
    rotor['pos'] = np.array([sval*rotor_dist, cval*rotor_dist, gearG['pos'][2]])

    for gear_name,gear in platter['gears'].items():
        if 'roll_to' in gear:
            roll_theta = np.radians(gear['roll_to'][0])
            roll_parent = platter['gears'][gear['roll_to'][1]]
            theta = gear.get('rot', 0.0)
            pos = gear.get('pos', [0.0, 0.0, 0.0])
            roll_radius = length(pos - roll_parent['pos'])
            sval = np.sin(roll_theta)
            cval = np.cos(roll_theta)
            gear['pos'] = np.array([sval*roll_radius + roll_parent['pos'][0], cval*roll_radius + roll_parent['pos'][1], gear['pos'][2]])
            gear['rot'] = theta + 2.0 * roll_theta * float(roll_parent['teeth']) / float(gear['teeth'])


    # set_platter_base(platter_name)
    # make_platter_supports(platter_name)
    # do_platter_adjustments(platter_name)

def print_checks():
    print('Pitch stats:')
    for platter_name,platter in all_platters.items():
        rotor_pin_width = None
        rotor = platter['gears'].get('Grotor', None)
        if rotor:
            rotor_pin_width = 2 * rotor['pin_radius']
        smallest_name = None
        smallest_pitch = 1000
        for gear_name,gear in platter['gears'].items():
            pitch = gear['specs'].get('pitch_ref', 1000)
            if pitch < smallest_pitch:
                smallest_pitch = pitch
                smallest_name = gear_name
        if smallest_name:
            print('  {}: {} {} rotor-notch: {}'.format(platter_name, smallest_name, smallest_pitch, rotor_pin_width))

class ColladaModel:
    def __init__(self, file_name):
        self.file_name = file_name
        self.start()

    def start(self):
        self.nodes = []
        self.materials = []
        self.mesh = col.Collada()

    def finish(self):
        myscene = col.scene.Scene("myscene", self.nodes)
        self.mesh.scenes.append(myscene)
        self.mesh.scene = myscene
        self.mesh.write(self.file_name)

    def new_material(self, diffuse_color, specular_color):
        mat_index = len(self.materials)
        effect = col.material.Effect("effect0"+str(mat_index), [], "phong", diffuse=diffuse_color, specular=specular_color)
        mat = col.material.Material("material0"+str(mat_index), "mymaterial"+str(mat_index), effect)
        self.mesh.effects.append(effect)
        self.mesh.materials.append(mat)
        self.materials.append(mat)
        return mat

    def add_extruded_tristrip_quad_mesh(self, material, verts, verts_z, pos=None, rot=None, closed=True):
        num_strip_verts = len(verts)
        num_sections = num_strip_verts >> 1
        num_quads = num_sections << 2
        vert_floats = np.ndarray((num_sections,4,3), dtype=np.float64)
        norm_floats = np.ndarray((num_sections,4,3), dtype=np.float64)
        indices     = np.ndarray((num_sections,4,6,2), dtype=np.uint32)
        z1 = np.array([0.0, 0.0, verts_z[1]])
        z2 = np.array([0.0, 0.0, verts_z[0]])
        for si in range(num_sections):
            vi = si << 1
            vert_floats[si][0] = verts[vi + 0] + z1
            vert_floats[si][1] = verts[vi + 0] + z2
            vert_floats[si][2] = verts[vi + 1] + z1
            vert_floats[si][3] = verts[vi + 1] + z2

            norm_floats[si][0] = np.array([0.0, 0.0,  1.0])
            norm_floats[si][1] = np.array([0.0, 0.0, -1.0])
            norm_floats[si][2] = np.array([0.0, 0.0, -1.0])
            norm_floats[si][3] = np.array([0.0, 0.0,  1.0])

            s0_start = 4 * si
            s1_start = 4 * ((si + 1) % num_sections)
            # top
            indices[si][0][0][0] = s0_start + 0
            indices[si][0][1][0] = s0_start + 2
            indices[si][0][2][0] = s1_start + 2
            indices[si][0][3][0] = s1_start + 2
            indices[si][0][4][0] = s1_start + 0
            indices[si][0][5][0] = s0_start + 0
            # bottom
            indices[si][1][0][0] = s0_start + 1
            indices[si][1][1][0] = s0_start + 3
            indices[si][1][2][0] = s1_start + 3
            indices[si][1][3][0] = s1_start + 3
            indices[si][1][4][0] = s1_start + 1
            indices[si][1][5][0] = s0_start + 1
            # in
            indices[si][2][0][0] = s0_start + 0
            indices[si][2][1][0] = s0_start + 1
            indices[si][2][2][0] = s1_start + 1
            indices[si][2][3][0] = s1_start + 1
            indices[si][2][4][0] = s1_start + 0
            indices[si][2][5][0] = s0_start + 0
            # out
            indices[si][3][0][0] = s0_start + 2
            indices[si][3][1][0] = s0_start + 3
            indices[si][3][2][0] = s1_start + 3
            indices[si][3][3][0] = s1_start + 3
            indices[si][3][4][0] = s1_start + 2
            indices[si][3][5][0] = s0_start + 2
            for section_face in range(4):
                ni = s0_start + section_face
                for nvi in range(6):
                    indices[si][section_face][nvi][1] = ni
    
        flat_verts   = vert_floats.flatten()
        flat_norms   = norm_floats.flatten()
        flat_indices = indices.flatten()

        node_index = len(self.nodes)
        vert_src = col.source.FloatSource("cubeverts-array"+str(node_index), np.array(flat_verts), ('X', 'Y', 'Z'))
        normal_src = col.source.FloatSource("cubenormals-array"+str(node_index), np.array(flat_norms), ('X', 'Y', 'Z'))
        geom = col.geometry.Geometry(self.mesh, "geometry0"+str(node_index), "mycube"+str(node_index), [vert_src, normal_src])

        input_list = col.source.InputList()
        input_list.addInput(0, 'VERTEX', "#cubeverts-array"+str(node_index))
        input_list.addInput(1, 'NORMAL', "#cubenormals-array"+str(node_index))

        triset = geom.createTriangleSet(indices, input_list, "materialref"+str(node_index))
        geom.primitives.append(triset)
        self.mesh.geometries.append(geom)

        matnode = col.scene.MaterialNode("materialref"+str(node_index), material, inputs=[])
        geomnode = col.scene.GeometryNode(geom, [matnode])
        node = col.scene.Node("node0"+str(node_index), children=[geomnode])
        self.nodes.append(node)


    def add_test_object(self, object_name, object_parent_name):
        vert_floats = [-50,50,50,50,50,50,-50,-50,50,50,
                       -50,50,-50,50,-50,50,50,-50,-50,-50,-50,50,-50,-50]
        normal_floats = [0,0,1,0,0,1,0,0,1,0,0,1,0,1,0,
                         0,1,0,0,1,0,0,1,0,0,-1,0,0,-1,0,0,-1,0,0,-1,0,-1,0,0,
                        -1,0,0,-1,0,0,-1,0,0,1,0,0,1,0,0,1,0,0,1,0,0,0,0,-1,
                        0,0,-1,0,0,-1,0,0,-1]
        vert_src = col.source.FloatSource("cubeverts-array", np.array(vert_floats), ('X', 'Y', 'Z'))
        normal_src = col.source.FloatSource("cubenormals-array", np.array(normal_floats), ('X', 'Y', 'Z'))
        geom = col.geometry.Geometry(self.mesh, "geometry0", "mycube", [vert_src, normal_src])

        input_list = col.source.InputList()
        input_list.addInput(0, 'VERTEX', "#cubeverts-array")
        input_list.addInput(1, 'NORMAL', "#cubenormals-array")

        indices = np.array([0,0,2,1,3,2,0,0,3,2,1,3,0,4,1,5,5,6,0,
                               4,5,6,4,7,6,8,7,9,3,10,6,8,3,10,2,11,0,12,
                               4,13,6,14,0,12,6,14,2,15,3,16,7,17,5,18,3,
                               16,5,18,1,19,5,20,7,21,6,22,5,20,6,22,4,23])
        triset = geom.createTriangleSet(indices, input_list, "materialref")
        geom.primitives.append(triset)
        self.mesh.geometries.append(geom)

        matnode = col.scene.MaterialNode("materialref", self.mat, inputs=[])
        geomnode = col.scene.GeometryNode(geom, [matnode])
        node = col.scene.Node("node0", children=[geomnode])
        self.nodes.append(node)



def build_all():
    platters_to_make = all_platters.keys()
    for platter_name in platters_to_make:
        make_platter_specs(platter_name)
    for platter_name in platters_to_make:
        build_one_platter(platter_name)
    for platter_name in platters_to_make:
        set_platter_base(platter_name)
    for platter_name in platters_to_make:
        make_platter_supports(platter_name)
    for platter_name in platters_to_make:
        do_platter_adjustments(platter_name)

    if use_collada:
        collada_model = ColladaModel('output/travelers_pocketwatch.dae')
    else:
        collada_model = None
    #collada_model.add_test_object('cube', 'cube_parent')

    for platter_name in platters_to_make:
        print('  Building platter {} ---------'.format(platter_name))
        write_stl_platter(platter_name, collada_model,
                          './output/platter_{}_drive_gears.stl'.format(platter_name),
                          './output/platter_{}_planet_gears.stl'.format(platter_name),
                          './output/platter_{}_frame.stl'.format(platter_name))
    if collada_model is not None:
        collada_model.finish()
    print_checks()

if __name__ == "__main__":
    build_all()
