import numpy as np
from ej_outer_tooth_shapes import outer_tooth_shapes_p10

# Global scale values
# All meassurements are in mm
spur_tooth_pitch = 2.0
geneva_tooth_pitch = 5.0
tooth_rim_thickness = 0.4
spur_thickness = 1.0
geneva_thickness = 1.0
platter_z = 4.0
g_rotor_pin_radius = 0.5
g_rotor_arm_length = 4.0
g_rotor_hub_radius = 3.0

"""
      best_errors [0.0002094545855868546, 0.002124126164289919, 0.0020335698811777547, 0.00014134032982582312]
      best_ratios [231.75420417869884, 1167.8433179723504, 1559.8970466919648, 59.061319340329824]
      best_gears M[(-27, 29, 13, 29, 23, -13, 21, 39),
                 V (-27, 29, 13, 31, 38, -29, 14, 38), 
                 M (-27, 29, 13, 11, 26, -28, 19, 37), 
                 L (-27, 29, 13, 23, 13, -7, 30, 37)]
      error_degrees_per_year [0.2376797302003289, 0.4783273624581241, -0.34284199857996345, 0.6293529646316143]
"""

all_platters = {
    'Drive': {
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':27},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':13},
            'Pr'     :{'type':'spur',   'inout':'out', 'teeth':29, 'outer_ring':1},
            'Pp'     :{'type':'spur',   'inout':'out', 'teeth':(29-13)>>1},
            'Grotor' :{'type':'rotor'},
        },
    },
    'Mercury': {
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':13},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':21},
            'Pr'     :{'type':'spur',   'inout':'out', 'teeth':39, 'outer_ring':1},
            'Pp'     :{'type':'spur',   'inout':'out', 'teeth':(39-21)>>1},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':29},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':23},
            'Grotor' :{'type':'rotor'},
        },
    },
    'Venus': {
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':29},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':14},
            'Pr'     :{'type':'spur',   'inout':'out', 'teeth':38, 'outer_ring':1},
            'Pp'     :{'type':'spur',   'inout':'out', 'teeth':(38-14)>>1},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':31},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':38},
            'Grotor' :{'type':'rotor'},
        },
    },
    'Mars': {
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':28},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':19},
            'Pr'     :{'type':'spur',   'inout':'out', 'teeth':37, 'outer_ring':1},
            'Pp'     :{'type':'spur',   'inout':'out', 'teeth':(37-19)>>1},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':11},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':26},
            'Grotor' :{'type':'rotor'},
        },
    },
    'Luna': {
        'gears': {
            'G'      :{'type':'geneva', 'inout':'out', 'teeth':7},
            'Ps'     :{'type':'spur',   'inout':'out', 'teeth':30},
            'Pp'     :{'type':'spur',   'inout':'out', 'teeth':37},
            'Da'     :{'type':'spur',   'inout':'out', 'teeth':23},
            'Db'     :{'type':'spur',   'inout':'out', 'teeth':13},
            'Grotor' :{'type':'rotor'},
        },
    },
}

def get_outer_tooth_shape(num_teeth, pitch_ref):
    original_shape = outer_tooth_shapes_p10[str(num_teeth)]
    shape = [0.1 * pitch_ref * np.array([-vert[0], vert[1], 0.0]) for vert in original_shape]
    return shape

def make_tooth_table(gear):
    # https://en.wikipedia.org/wiki/Involute_gear
    # Symbols from https://khkgears.net/new/gear_knowledge/gear_technical_reference/involute_gear_profile.html
    num_teeth = gear['teeth']
    pitch_ref = spur_tooth_pitch

    p = pitch_ref
    m = p / np.pi # Module is the unit size indicated in millimeter (mm)
    ha = 1.0 * m # The distance between reference line and tooth tip
    hf = 1.25 * m # The distance between reference line and tooth root
    h = ha + hf # The distance between tooth tip and tooth root
    hw = 2.0 * m # Depth of tooth meshed with the mating gear
    c = 0.25 * m # The distance (clearance) between tooth root and the tooth tip of mating gear
    rhof = 0.38 * m # The radius of curvature between tooth surface and the tooth root
    rc = p * num_teeth / (2.0 * np.pi) # circular (reference) radius
    ri = rc - hf # inner (root) radius
    ro = rc + ha # outer (tip) radius

    gear['specs'] = {
                     'pitch_ref':pitch_ref,
                     'radius_ref':rc,
                     'radius_inner':ri,
                     'radius_outer':ro,
                     }
    face_segments = 10
    pitch_theta = 2.0 * np.pi / num_teeth
    print('num_teeth', num_teeth)
    print('radius_ref', rc, 'radius_inner', ri, 'radius_outer', ro)
    print('pitch_theta', pitch_theta)
    print('')
    in_theta = 0.0
    ref_theta = get_involute_theta(rc, ri)
    out_theta = get_involute_theta(ro, ri)
    total_tooth_theta = 2.0 * np.pi / num_teeth
    inner_flat_theta = 0.5 * total_tooth_theta - 2.0 * ref_theta
    outer_flat_theta = 0.5 * total_tooth_theta - 2.0 * (out_theta - ref_theta)
    print('ref_theta',ref_theta)
    print('out_theta',out_theta)
    print('inner_flat_theta',inner_flat_theta)
    print('outer_flat_theta',outer_flat_theta)
    print('total_tooth_theta',total_tooth_theta)

    tooth_curve_segments = 10
    # tc_radius1 = []
    # tc_theta1  = []
    # for i in range(tooth_curve_segments):
    #     t = float(i) / tooth_curve_segments
    #     r = (rc - ri) * t + ri
    #     theta = get_involute_theta(r, ri)
    #     tc_radius1.append(r)
    #     tc_theta1.append(theta)
    # tc_radius2 = []
    # tc_theta2  = []
    # for i in range(tooth_curve_segments):
    #     t = float(i) / tooth_curve_segments
    #     r = (ro - rc) * t + rc
    #     theta = get_involute_theta(r, ri)
    #     tc_radius2.append(r)
    #     tc_theta2.append(theta)
    result_theta = []
    result_radius = []
    # Inner flat
    theta_start = 0.0
    for i in range(tooth_curve_segments):
        t = float(i) / tooth_curve_segments
        theta = inner_flat_theta * t
        result_theta.append(theta + theta_start)
        result_radius.append(ri)
    theta_start += inner_flat_theta
    # Tooth rise
    for i in range(tooth_curve_segments):
        t = float(i) / tooth_curve_segments
        r = (ro - ri) * t + ri
        theta = get_involute_theta(r, ri)
        result_theta.append(theta + theta_start)
        result_radius.append(r)
    theta_start += out_theta
    # Outer flat
    for i in range(tooth_curve_segments):
        t = float(i) / tooth_curve_segments
        theta = outer_flat_theta * t
        result_theta.append(theta + theta_start)
        result_radius.append(ro)
    theta_start += outer_flat_theta
    # Tooth fall
    for i in range(tooth_curve_segments):
        t = float(i) / tooth_curve_segments
        r = (ro - ri) * (1.0 - t) + ri
        theta = get_involute_theta(r, ri)
        result_theta.append((out_theta - theta) + theta_start)
        result_radius.append(r)


    return[result_theta, result_radius]

def get_involute_theta(radius, root_radius):
    alpha = np.arccos(root_radius / radius)
    involute_alpha = np.tan(alpha) - alpha
    return np.arcsin(involute_alpha)

####### This is *almost* right, but
####### the correect version can be found here: http://hessmer.org/gears/InvoluteSpurGearBuilder.html
def build_one_spur(gear):
    num_teeth = gear['teeth']
    pitch_ref = spur_tooth_pitch
    gear['specs'] = {
                     'pitch_ref':pitch_ref,
                     'radius_ref':pitch_ref * num_teeth / (2.0 * np.pi),
                     # 'radius_inner':ri,
                     # 'radius_outer':ro,
                     }
    specs = gear['specs']
    tooth_table2 = get_outer_tooth_shape(num_teeth, specs['pitch_ref'])
    num_verts = 2 * num_teeth * len(tooth_table2)
    verts = np.ndarray((num_verts, 3), dtype=np.float64)

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
            outer_vert_index = tooth_vert_index + i*2
            verts[outer_vert_index][0] = tx * cval + ty * sval
            verts[outer_vert_index][1] = ty * cval - tx * sval
            verts[outer_vert_index][2] = 0.0

    outer_mult = -1.0 if gear.get('outer_ring', 0) else 1.0
    # Inner verts
    for tooth in range(num_teeth):
        tooth_vert_index = 2 * tooth * len(tooth_table2)
        for i in range(len(tooth_table2)):
            outer_vert_index = tooth_vert_index + i*2
            inner_vert_index = outer_vert_index + 1
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
    return verts

def make_geneva_tooth_table(gear):
    num_teeth = gear['teeth']
    pitch_ref = geneva_tooth_pitch
    gear['specs'] = {
                     'pitch_ref':pitch_ref,
                     'radius_ref':pitch_ref * num_teeth / (2.0 * np.pi),
                     'radius_inner':1.0,
                     'radius_outer':2.0,
                     }
    r = gear['specs']['radius_ref']
    total_tooth_theta = 2.0 * np.pi / num_teeth
    tooth_segments = 20
    result_theta = [total_tooth_theta * float(x) / float(tooth_segments) for x in range(tooth_segments)]
    result_radius = [r]*tooth_segments
    return [result_theta, result_radius]

def build_one_geneva(gear):
    tooth_table = make_geneva_tooth_table(gear)
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
            d1 = prev_outer_vert - this_outer_vert
            d2 = next_outer_vert - this_outer_vert
            c = np.cross(d1, d2)
            inner_vert_dir = np.sign(c[2]) * normalized(d1 + d2)
            verts[inner_vert_index] = this_outer_vert + tooth_rim_thickness * inner_vert_dir
    gear['verts'] = [verts]
    return verts

def build_one_rotor(gear):
    g_rotor_hub_inner_radius = 0.1 * g_rotor_hub_radius
    num_segments = 400
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
        outer_radius = g_rotor_hub_radius
        inner_radius = g_rotor_hub_inner_radius
        outer_vert_index = seg*2
        inner_vert_index = outer_vert_index + 1
        verts[outer_vert_index][0] = outer_radius * sval
        verts[outer_vert_index][1] = outer_radius * cval
        verts[outer_vert_index][2] = 0.0
        verts[inner_vert_index][0] = inner_radius * sval
        verts[inner_vert_index][1] = inner_radius * cval
        verts[inner_vert_index][2] = 0.0

        outer_pin_radius = g_rotor_pin_radius
        inner_pin_radius = g_rotor_pin_radius * 0.1
        pin_verts[outer_vert_index][0] = g_rotor_arm_length + outer_pin_radius * sval
        pin_verts[outer_vert_index][1] = outer_pin_radius * cval
        pin_verts[outer_vert_index][2] = 0.0
        pin_verts[inner_vert_index][0] = g_rotor_arm_length + inner_pin_radius * sval
        pin_verts[inner_vert_index][1] = inner_pin_radius * cval
        pin_verts[inner_vert_index][2] = 0.0

        spur_thickness = 1.0
        outer_disc_radius = 2.0 * g_rotor_pin_radius + g_rotor_arm_length
        inner_disc_radius = g_rotor_hub_inner_radius
        disc_verts[outer_vert_index][0] = outer_disc_radius * sval
        disc_verts[outer_vert_index][1] = outer_disc_radius * cval
        disc_verts[outer_vert_index][2] = -1.0 * spur_thickness
        disc_verts[inner_vert_index][0] = inner_disc_radius * sval
        disc_verts[inner_vert_index][1] = inner_disc_radius * cval
        disc_verts[inner_vert_index][2] = -1.0 * spur_thickness

    gear['verts'] = [verts, pin_verts, disc_verts]
    return verts

def normalized(v):
    return v / np.linalg.norm(v)

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

def write_stl_tristrip_quads(fp, verts, pos=None, rot=None, closed=True):
    if pos is None:
        pos = np.array([0,0,0])
    if rot is None:
        rot = 0.0
    sval = np.sin(rot)
    cval = np.cos(rot)
    xverts = verts + pos
    num_verts = len(verts)
    num_quads = num_verts >> 1
    if not closed:
        num_quads -= 1
    spur_thickness = 1.0
    gear_half_z = np.array([0.0, 0.0, 0.5 * spur_thickness])
    n_up = np.array([0.0, 0.0, 1])
    n_down = np.array([0.0, 0.0, -1])
    for quad in range(num_quads):
        v0a = xverts[quad * 2] + gear_half_z
        v1a = xverts[quad * 2 + 1] + gear_half_z
        v2a = xverts[(quad * 2 + 2) % num_verts] + gear_half_z
        v3a = xverts[(quad * 2 + 3) % num_verts] + gear_half_z
        v0b = v0a - 2 * gear_half_z
        v1b = v1a - 2 * gear_half_z
        v2b = v2a - 2 * gear_half_z
        v3b = v3a - 2 * gear_half_z
        write_one_quad(fp, v0a, v1a, v2a, v3a, n=n_up)
        write_one_quad(fp, v1b, v0b, v3b, v2b, n=n_down)
        write_one_quad(fp, v0b, v0a, v2b, v2a)
        write_one_quad(fp, v1a, v1b, v3a, v3b)

def write_stl_platter(platter_name, file_name):
    platter = all_platters[platter_name]
    # try:
    fp = open(file_name, 'w')
    for gear_name,gear in platter['gears'].items():
        print('  Writing gear {}:{}'.format(platter_name, gear_name))
        fp.write('solid OpenSCAD_Model\n')
        verts = gear.get('verts', None)
        if verts is not None:
            for strip_verts in verts:
                pos = gear.get('pos', None)
                rot = gear.get('rot', None)
                write_stl_tristrip_quads(fp, strip_verts, pos=pos, rot=rot, closed=True)
        fp.write('endsolid OpenSCAD_Model\n')
    fp.close()
    # except:
    #     print('Unable to write file', file_name)

def build_one_platter(platter_name):
    platter = all_platters[platter_name]
    gearG  = platter['gears']['G']
    gearPr = platter['gears']['Pr']
    gearPp = platter['gears']['Pp']
    gearPs = platter['gears']['Ps']
    gearDa = platter['gears']['Da']
    gearDb = platter['gears']['Db']
    rotor = platter['gears']['Grotor']

    print('  Building platter {}'.format(platter_name))
    for gear_name,gear in platter['gears'].items():
        gtype = gear['type']
        if gtype == 'geneva':
            build_one_geneva(gear)
        elif gtype == 'spur':
            build_one_spur(gear)
        elif gtype == 'rotor':
            build_one_rotor(gear)
    gearG['pos'] = np.array([0.0, 0.0, 0.0])
    gearPs['pos'] = np.array([0.0, 0.0, platter_z])
    gearPr['pos'] = np.array([0.0, 0.0, platter_z])
    gearPp['pos'] = np.array([gearPs['specs']['radius_ref'] + gearPp['specs']['radius_ref'], 0.0, platter_z])
    gearDa['pos'] = np.array([0.0, 0.0, -platter_z])
    gearDb['pos'] = np.array([0.0, gearDa['specs']['radius_ref'] + gearDb['specs']['radius_ref'], -platter_z])
    rotor['pos'] = np.array([gearDb['pos'][0], gearDb['pos'][1], gearG['pos'][2]])

def build_all():

    build_one_platter('Venus')
    write_stl_platter('Venus', './output/platter_venus.stl')

if __name__ == "__main__":
    build_all()
