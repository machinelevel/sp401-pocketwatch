###

km_per_au = 149598000

mercury = {'name':'mercury', 'nspace':'', 'orbital_period_days': 87.9691,       'sydonic_period_days':0,          'min_dist_au': 0.307499, 'max_dist_au': 0.466697}
venus   = {'name':'venus', 'nspace':'  ', 'orbital_period_days': 224.701,       'sydonic_period_days':0,          'min_dist_au': 0.718440, 'max_dist_au': 0.728213}
earth   = {'name':'earth', 'nspace':'  ', 'orbital_period_days': 365.256363004, 'sydonic_period_days':0,          'min_dist_au': 0.98327,  'max_dist_au': 1.017}
mars    = {'name':'mars', 'nspace':'   ', 'orbital_period_days': 686.971,       'sydonic_period_days':0,          'min_dist_au': 1.382,    'max_dist_au': 1.666}
jupiter = {'name':'jupiter', 'nspace':'', 'orbital_period_days': 0,             'sydonic_period_days':0,          'min_dist_au': 0,        'max_dist_au': 0}
saturn  = {'name':'saturn', 'nspace':' ', 'orbital_period_days': 0,             'sydonic_period_days':0,          'min_dist_au': 0,        'max_dist_au': 0}
uranus  = {'name':'uranus', 'nspace':' ', 'orbital_period_days': 0,             'sydonic_period_days':0,          'min_dist_au': 0,        'max_dist_au': 0}
neptune = {'name':'neptune', 'nspace':'', 'orbital_period_days': 0,             'sydonic_period_days':0,          'min_dist_au': 0,        'max_dist_au': 0}
pluto   = {'name':'pluto',  'nspace':' ', 'orbital_period_days': 0,             'sydonic_period_days':0,          'min_dist_au': 0,        'max_dist_au': 0}
luna    = {'name':'luna', 'nspace':'   ', 'orbital_period_days': 27.321661,     'sydonic_period_days':29.530589,  'min_dist_au': 0.00257,              'max_dist_au': 0.00257}
phobos  = {'name':'phobos', 'nspace':' ', 'orbital_period_days': 0.31891023,    'sydonic_period_days':0.31891023, 'min_dist_au': 9376 / km_per_au,     'max_dist_au': 9376 / km_per_au}
deimos  = {'name':'deimos', 'nspace':' ', 'orbital_period_days': 1.263,         'sydonic_period_days':1.263,      'min_dist_au': 23463.2 / km_per_au,  'max_dist_au': 23463.2 / km_per_au}

planets = [mercury, venus, earth, mars]#, jupiter, saturn, uranus, neptune, pluto]
planets_without_earth = [mercury, venus, mars]#, jupiter, saturn, uranus, neptune, pluto]
moons =   [luna]

def list_planets():
    for planet in planets:
        print('  {}{}: orbital period {} days, {} orbits per earth orbit'.format(planet['nspace'], planet['name'],
                 planet['orbital_period_days'],
                 earth['orbital_period_days'] / planet['orbital_period_days']))

def make_conjunctions():
    print('conjunctions ---------')
    for planet in planets:
        if planet is not earth:
            orbit_per_earth_year = earth['orbital_period_days'] / planet['orbital_period_days']
            relative_orbit_per_earth_year = orbit_per_earth_year - 1.0
            earth_years_per_relative_orbit = 1.0 / relative_orbit_per_earth_year
            days_per_relative_orbit = earth_years_per_relative_orbit * earth['orbital_period_days']
            planet['sydonic_period_days'] = days_per_relative_orbit
            print('  {}{}: sydonic_period_days {}'.format(planet['nspace'], planet['name'], planet['sydonic_period_days']))
    for moon in moons:
        print('  {}{}: sydonic_period_days {}'.format(moon['nspace'], moon['name'], moon['sydonic_period_days']))

def make_gear_ratios():
    print('mobile earth ---------')
    for planet in planets:
        planet['em_gear_ratio'] = 2.0 * planet['orbital_period_days']
        print('  {}{}: hour-hand-ratio 1:{}'.format(planet['nspace'], planet['name'], planet['em_gear_ratio']))
    print('static earth ---------')
    for planet in planets:
        if planet is not earth:
            planet['es_gear_ratio'] = 2.0 * planet['sydonic_period_days']
            print('  {}{}: hour-hand-ratio 1:{}'.format(planet['nspace'], planet['name'], planet['es_gear_ratio']))
    for moon in moons:
        moon['es_gear_ratio'] = 2.0 * moon['sydonic_period_days']
        print('  {}{}: hour-hand-ratio 1:{}'.format(moon['nspace'], moon['name'], moon['es_gear_ratio']))

def print_gear_ratio_search_progress(bodies, gears, best_errors, best_ratios, best_gears, is_sydonic):
    print('  searching {}'.format(gears))
    print('      best_errors {}'.format(best_errors))
    print('      best_ratios {}'.format(best_ratios))
    print('      best_gears {}'.format(best_gears))
    error_degrees_per_year = [100]*len(best_errors)
    if is_sydonic:
        for i,body in enumerate(bodies):
            error = best_errors[i]
            days_per_relative_orbit = body['sydonic_period_days']
            earth_years_per_relative_orbit = days_per_relative_orbit / earth['orbital_period_days']
            relative_orbit_per_earth_year = 1.0 / earth_years_per_relative_orbit
            error_degrees_per_year[i] = error * relative_orbit_per_earth_year * 360.0
    else:
        for i,body in enumerate(bodies):
            error = best_errors[i]
            days_per_orbit = body['orbital_period_days']
            earth_years_per_orbit = days_per_orbit / earth['orbital_period_days']
            orbit_per_earth_year = 1.0 / earth_years_per_orbit
            error_degrees_per_year[i] = error * orbit_per_earth_year * 360.0

    print('      error_degrees_per_year {}'.format(error_degrees_per_year))

def old_gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic):
    # if chain_name == 'chain-v2':
    #     gear_ratio_search_chain_v2(chain_name, bodies, desired_ratios, is_sydonic)
    g1_range = range(1, 30+1)
    g2_range = range(1, 30+1)
    c1_range = range(8, 40+1)
    c2_range = range(8, 40+1)
    a_range = range(8, 40+1)
    s_range = range(8, 40+1)
    g1_range = [19]
#    g2_range = [22]
    best_errors = [1.0] * len(desired_ratios)
    best_ratios = [1.0] * len(desired_ratios)
    best_gears = [(1,1,1,1)] * len(desired_ratios)
    best_drive = ['(none)'] * len(desired_ratios)
    do_planetaries = False
    for g1 in g1_range:
        for g2 in g2_range:
            if 1 or g2 <= g1:
                print_gear_ratio_search_progress(bodies, (g1, g2), best_errors, best_ratios, best_gears, is_sydonic)
                for c1 in c1_range:
                    for c2 in c2_range:
                        for a in a_range:
                            for s in s_range:
                                correction_drive = 'a->s'
                                correction_ratio = a / s
                                ratio = g1 * g2 * (c2 / c1) * correction_ratio
                                for i,dr in enumerate(desired_ratios):
                                    error = abs(abs(dr) - ratio)
                                    if error < best_errors[i]:
                                        best_errors[i] = error
                                        best_gears[i] = (g1, g2, c1, c2, a, s)
                                        best_ratios[i] = ratio
                                        best_drive[i] = correction_drive
                                if a > s and do_planetaries:
                                    correction_drive = 'a->p'
                                    correction_ratio = 1 + (a / s)
                                    ratio = g1 * g2 * (c2 / c1) * correction_ratio
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (g1, g2, c1, c2, a, s)
                                            best_ratios[i] = ratio
                                            best_drive[i] = correction_drive
                                    correction_drive = 'p->a'
                                    correction_ratio = (s + a) / a
                                    ratio = g1 * g2 * (c2 / c1) * correction_ratio
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (g1, g2, c1, c2, a, s)
                                            best_ratios[i] = ratio
                                            best_drive[i] = correction_drive
                                    correction_drive = 'p->s'
                                    correction_ratio = 1.0 / (1 + (a / s))
                                    ratio = g1 * g2 * (c2 / c1) * correction_ratio
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (g1, g2, c1, c2, a, s)
                                            best_ratios[i] = ratio
                                            best_drive[i] = correction_drive
                                    correction_drive = 's->p'
                                    correction_ratio = 1.0 / ((s + a) / a)
                                    ratio = g1 * g2 * (c2 / c1) * correction_ratio
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (g1, g2, c1, c2, a, s)
                                            best_ratios[i] = ratio
                                            best_drive[i] = correction_drive

    return({'desired_ratios':desired_ratios, 'best_drive':best_drive, 'best_errors':best_errors, 'best_ratios':best_ratios, 'best_gears':best_gears})

def gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic):
# Here are the first best v2 numbers
# desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
# best_drive [None, None, None, None]
# best_errors [0.0003813752198880138, 0.0002808081917464733, 0.0025200847474025068, 0.0005503950617153919]
# best_ratios [231.75403225806454, 1167.845161290323, 1559.897533206831, 59.061728395061714]
# best_gears [(10, 20, 31, 11, -19, -25, 11, 16), (10, 20, 10, 22, -17, -22, 22, 31), (10, 20, 17, 37, -21, -23, 23, 31), (10, 20, 27, 13, -23, 9, 24)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [0.43276760503525186, 0.06323458746134669, -0.42486412657231726, 2.450770875064815]

# This is also good (except merc) considering the number of shared gears:
# desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
# best_drive [None, None, None, None]
# best_errors [0.05070992958079046, 0.0339942870671166, 0.03257925964999231, 0.0005503950617153919]
# best_ratios [231.70370370370364, 1167.8114478114476, 1559.8624338624336, 59.061728395061714]
# best_gears [(10, 20, 27, 13, -23, -24, 17, 39), (10, 20, 27, 13, -23, -29, 20, 11), (10, 20, 27, 13, -23, -29, 34, 14), (10, 20, 27, 13, -23, 9, 24)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [57.54336839878852, 7.655099751051719, -5.49257667220662, 2.450770875064815]

# This is also good and V and M are such slight corrections maybe they're not needed?
# desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
# best_drive [None, None, None, None]
# best_errors [0.0006009971188234431, 0.0013875176180135895, 0.018469912207137895, 0.0005503950617153919]
# best_ratios [231.7538126361656, 1167.8440545808967, 1559.8765432098764, 59.061728395061714]
# best_gears [(10, 20, 27, 23, -5, -25, 37, 34), (10, 20, 27, 32, -22, -23, 37, 38), (10, 20, 27, 38, -19, -30, 35, 36), (10, 20, 27, 13, -23, 9, 24)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [0.6819847493570269, 0.31245208205911135, -3.1138647721435815, 2.450770875064815]

    # if chain_name == 'chain-v2':
    #     gear_ratio_search_chain_v2(chain_name, bodies, desired_ratios, is_sydonic)
    gHr_range = [10]
    g24a_range = [20] # This is the 24h rotator
    g24b_range = range(8, 40+1) # drives G1D
    gGDa_range = range(8, 40+1) # Driver to Geneva1
    gG1_range = range(1, 30+1) # The first real Geneva
    gG2_range = range(1, 30+1) # The second real Geneva
    gG1L_range = range(8, 40+1) # G1 to Luna
    gL_range = range(8, 40+1) # Luna
    # gG2D_range = range(8, 40+1)
    gG2DMerc_range = range(8, 40+1)
    gG2DVenus_range = range(8, 40+1)
    gG2DMars_range = range(8, 40+1)
    gMerc_range = range(8, 40+1)
    gVenus_range = range(8, 40+1)
    gMars_range = range(8, 40+1)

    # Constraints for gear sharing
    g24b_range = [27]
    gGDa_range = [38]
    gG1_range = [19]
#    gG2_range = [23]

#    constraints: 
#   g24a = 2 * gHr
#   gG2D all the same
#   gG2 >= gG1 (since they're interchangeable)

    best_errors = [1.0] * len(desired_ratios)
    best_ratios = [1.0] * len(desired_ratios)
    best_gears = [None] * len(desired_ratios)
    best_drive = ['(none)'] * len(desired_ratios)
    do_planetaries = False
    iLuna = len(bodies) - 1
    iMerc = 0
    iVenus = 1
    iMars = 2
    for gHr in gHr_range:
        g24_ratio = 2.0
        g24a = 2 * gHr
        for g24b in g24b_range:
            for gGDa in gGDa_range:
                gGD_ratio = g24_ratio * gGDa / g24b
                print_gear_ratio_search_progress(bodies, (g24a, g24b, gGDa), best_errors, best_ratios, best_gears, is_sydonic)
                for gG1 in gG1_range:
                    gG1_ratio = gG1 * gGD_ratio

                    for gG1L in gG1L_range:
                        for gL in gL_range:
                            gL_ratio = gG1_ratio * gL / gG1L
                            i = iLuna
                            dr = desired_ratios[i]
                            ratio = gL_ratio
                            error = abs(abs(dr) - ratio)
                            if error < best_errors[i]:
                                best_errors[i] = error
                                best_gears[i] = (gHr, g24a, g24b, gGDa, -gG1, gG1L, gL)
                                best_ratios[i] = ratio
                                best_drive[i] = None

                    for gG2 in gG2_range:
                        if gG2 >= gG1:
                            gG2_ratio = gG1_ratio * gG2
                            for gG2D in gG2DMars_range:
                                gG2D_ratio = gG2_ratio * gG2D
                                for gMars in gMars_range:
                                    if gMars > gG2D:
                                        gMars_ratio = gG2D_ratio / gMars
                                        i = iMars
                                        dr = desired_ratios[i]
                                        ratio = gMars_ratio
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (gHr, g24a, g24b, gGDa, -gG1, -gG2, gG2D, gMars)
                                            best_ratios[i] = ratio
                                            best_drive[i] = None
                            for gG2D in gG2DMerc_range:
                                gG2D_ratio = gG2_ratio * gG2D
                                for gMerc in gMerc_range:
                                    gMerc_ratio = gG2D_ratio / gMerc
                                    i = iMerc
                                    dr = desired_ratios[i]
                                    ratio = gMerc_ratio
                                    error = abs(abs(dr) - ratio)
                                    if error < best_errors[i]:
                                        best_errors[i] = error
                                        best_gears[i] = (gHr, g24a, g24b, gGDa, -gG1, -gG2, gG2D, gMerc)
                                        best_ratios[i] = ratio
                                        best_drive[i] = None
                            for gG2D in gG2DVenus_range:
                                gG2D_ratio = gG2_ratio * gG2D
                                for gVenus in gVenus_range:
                                    if gVenus > gG2D:
                                        gVenus_ratio = gG2D_ratio / gVenus
                                        i = iVenus
                                        dr = desired_ratios[i]
                                        ratio = gVenus_ratio
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (gHr, g24a, g24b, gGDa, -gG1, -gG2, gG2D, gVenus)
                                            best_ratios[i] = ratio
                                            best_drive[i] = None

    return({'desired_ratios':desired_ratios, 'best_drive':best_drive, 'best_errors':best_errors, 'best_ratios':best_ratios, 'best_gears':best_gears})


            # if 1 or g2 <= g1:
            #     print_gear_ratio_search_progress(bodies, (g1, g2), best_errors, best_ratios, best_gears, is_sydonic)
            #     for c1 in c1_range:
            #         for c2 in c2_range:
            #             for a in a_range:
            #                 for s in s_range:
            #                     correction_drive = 'a->s'
            #                     correction_ratio = a / s
            #                     ratio = g1 * g2 * (c2 / c1) * correction_ratio
            #                     for i,dr in enumerate(desired_ratios):
            #                         error = abs(abs(dr) - ratio)
            #                         if error < best_errors[i]:
            #                             best_errors[i] = error
            #                             best_gears[i] = (g1, g2, c1, c2, a, s)
            #                             best_ratios[i] = ratio
            #                             best_drive[i] = correction_drive
            #                     if a > s and do_planetaries:
            #                         correction_drive = 'a->p'
            #                         correction_ratio = 1 + (a / s)
            #                         ratio = g1 * g2 * (c2 / c1) * correction_ratio
            #                         for i,dr in enumerate(desired_ratios):
            #                             error = abs(abs(dr) - ratio)
            #                             if error < best_errors[i]:
            #                                 best_errors[i] = error
            #                                 best_gears[i] = (g1, g2, c1, c2, a, s)
            #                                 best_ratios[i] = ratio
            #                                 best_drive[i] = correction_drive
            #                         correction_drive = 'p->a'
            #                         correction_ratio = (s + a) / a
            #                         ratio = g1 * g2 * (c2 / c1) * correction_ratio
            #                         for i,dr in enumerate(desired_ratios):
            #                             error = abs(abs(dr) - ratio)
            #                             if error < best_errors[i]:
            #                                 best_errors[i] = error
            #                                 best_gears[i] = (g1, g2, c1, c2, a, s)
            #                                 best_ratios[i] = ratio
            #                                 best_drive[i] = correction_drive
            #                         correction_drive = 'p->s'
            #                         correction_ratio = 1.0 / (1 + (a / s))
            #                         ratio = g1 * g2 * (c2 / c1) * correction_ratio
            #                         for i,dr in enumerate(desired_ratios):
            #                             error = abs(abs(dr) - ratio)
            #                             if error < best_errors[i]:
            #                                 best_errors[i] = error
            #                                 best_gears[i] = (g1, g2, c1, c2, a, s)
            #                                 best_ratios[i] = ratio
            #                                 best_drive[i] = correction_drive
            #                         correction_drive = 's->p'
            #                         correction_ratio = 1.0 / ((s + a) / a)
            #                         ratio = g1 * g2 * (c2 / c1) * correction_ratio
            #                         for i,dr in enumerate(desired_ratios):
            #                             error = abs(abs(dr) - ratio)
            #                             if error < best_errors[i]:
            #                                 best_errors[i] = error
            #                                 best_gears[i] = (g1, g2, c1, c2, a, s)
            #                                 best_ratios[i] = ratio
            #                                 best_drive[i] = correction_drive


def snowman_gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic):
    """
    These are great!
    Here's the best with drive of 9:

    desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
    best_drive [None, None, None, None]
    best_errors [0.0013272135313400213, 0.012108765181437775, 0.010955151069083513, 0.0005503950617296027]
    best_ratios [231.7530864197531, 1167.8333333333333, 1559.8840579710145, 59.06172839506173]
    best_gears [(9, 8, -13, -19, 18, 19), (9, 11, -7, -7, 2, 39), (9, 31, -12, -28, 23, 31), (9, 8, -1, -23, 9, 26)]
    Venus G2 *= 4 so Venus G2b isn't 2...
    best_gears [(9, 8, -13, -19, 18, 19), (9, 11, -7, -28, 8, 39), (9, 31, -12, -28, 23, 31), (9, 8, -1, -23, 9, 26)]
    bodies ['mercury', '  venus', '   mars', '   luna']
    error_degrees_per_year [1.5060627732894087, 2.7267465601781287, -1.8469421297166195, 2.4507708751280926]

    ...but this is also very good, with smaller Mars gear:
      best_errors [0.0013272135313400213, 0.012108765181437775, 0.018469912207137895, 0.0005503950617296027]
      best_ratios [231.7530864197531, 1167.8333333333333, 1559.8765432098764, 59.06172839506173]
      best_gears [(9, 8, -13, -19, 18, 19), (9, 11, -7, -28, 8, 39), (9, 10, -19, -19, 9, 35), (9, 8, -1, -23, 9, 26)]
      error_degrees_per_year [1.5060627732894087, 2.7267465601781287, -3.1138647721435815, 2.4507708751280926]

best_gears [(9, 8, -13, -19, 18, 19), (9, 11, -7, -28, 8, 39), (9, 31, -12, -28, 23, 31), (9, 8, -1, -23, 9, 26)]
bodies ['mercury', '  venus', '   mars', '   luna']
error_degrees_per_year [1.5060627732894087, 2.7267465601781287, -1.8469421297166195, 2.4507708751280926]

    """
    gHr_range = range(8, 20+1)
    gMin_range = range(8, 20+1)
    g24a_range = range(16, 40+1) # This is the 24h rotator

    gMerc1_range  = range(8, 40+1) # Connects to Hr or Min
    gVenus1_range = range(8, 40+1) # Connects to Hr or Min
    gMars1_range  = range(8, 40+1) # Connects to Hr or Min
    gVenus1_range = range(8, 40+1) # Connects to Hr or Min
    gLuna1_range  = range(8, 40+1) # Connects to Hr or Min

    gMercG1_range   = range(1, 30+1) # Geneva stage 1
    gVenusG1_range  = range(1, 30+1) # Geneva stage 1
    gMarsG1_range   = range(1, 30+1) # Geneva stage 1
    gVenusG1_range  = range(1, 30+1) # Geneva stage 1
    gLunaG1_range   = range(1, 30+1) # Geneva stage 1

    gMercG2_range   = range(1, 30+1) # Geneva stage 2
    gVenusG2_range  = range(1, 30+1) # Geneva stage 2
    gMarsG2_range   = range(1, 30+1) # Geneva stage 2
    gVenusG2_range  = range(1, 30+1) # Geneva stage 2
    gLunaG2_range   = range(1, 30+1) # Geneva stage 2

    gMercG2b_range   = range(8, 40+1) # teeth on G2
    gVenusG2b_range  = range(8, 40+1) # teeth on G2
    gMarsG2b_range   = range(8, 40+1) # teeth on G2
    gVenusG2b_range  = range(8, 40+1) # teeth on G2
    gLunaG2b_range   = range(8, 40+1) # teeth on G2

    gMercD_range   = range(8, 40+1) # teeth on driver
    gVenusD_range  = range(8, 40+1) # teeth on driver
    gMarsD_range   = range(8, 40+1) # teeth on driver
    gVenusD_range  = range(8, 40+1) # teeth on driver
    gLunaD_range   = range(8, 40+1) # teeth on driver

    gHr_range = [18]


    best_errors = [1.0] * len(desired_ratios)
    best_ratios = [1.0] * len(desired_ratios)
    best_gears = [None] * len(desired_ratios)
    best_drive = ['(none)'] * len(desired_ratios)
    do_planetaries = False
    iLuna = len(bodies) - 1
    iMerc = 0
    iVenus = 1
    iMars = 2
    for gHr in gHr_range:
        g24_ratio = 2.0
        g24a = g24_ratio * gHr
        for gMerc1 in gMerc1_range:
            for gMercG1 in gMercG1_range:
                for gMercG2 in gMercG2_range:
                    print_gear_ratio_search_progress(bodies, (gHr, gMerc1, gMercG1, gMercG2), best_errors, best_ratios, best_gears, is_sydonic)
                    for gMercG2b in gMercG2b_range:
                        for gMercD in gMercD_range:
                            ratio = (gMerc1 * gMercG1 * gMercG2 * gMercD) / (gHr * gMercG2b)
                            for i,dr in enumerate(desired_ratios):
                                error = abs(abs(dr) - ratio)
                                if error < best_errors[i]:
                                    best_errors[i] = error
                                    best_gears[i] = (gHr, gMerc1, -gMercG1, -gMercG2, gMercG2b, gMercD)
                                    best_ratios[i] = ratio
                                    best_drive[i] = None
                            # @@@@@ find a simpler way to put planet-carrier in the search
                            # a = gMercG2b
                            # if a > s and do_planetaries:
                            #     correction_drive = 'a->p'
                            #     correction_ratio = 1 + (a / s)
                            #     ratio = g1 * g2 * (c2 / c1) * correction_ratio
                            #     for i,dr in enumerate(desired_ratios):
                            #         error = abs(abs(dr) - ratio)
                            #         if error < best_errors[i]:
                            #             best_errors[i] = error
                            #             best_gears[i] = (g1, g2, c1, c2, a, s)
                            #             best_ratios[i] = ratio
                            #             best_drive[i] = correction_drive
                            #     correction_drive = 'p->a'
                            #     correction_ratio = (s + a) / a
                            #     ratio = g1 * g2 * (c2 / c1) * correction_ratio
                            #     for i,dr in enumerate(desired_ratios):
                            #         error = abs(abs(dr) - ratio)
                            #         if error < best_errors[i]:
                            #             best_errors[i] = error
                            #             best_gears[i] = (g1, g2, c1, c2, a, s)
                            #             best_ratios[i] = ratio
                            #             best_drive[i] = correction_drive
                            #     correction_drive = 'p->s'
                            #     correction_ratio = 1.0 / (1 + (a / s))
                            #     ratio = g1 * g2 * (c2 / c1) * correction_ratio
                            #     for i,dr in enumerate(desired_ratios):
                            #         error = abs(abs(dr) - ratio)
                            #         if error < best_errors[i]:
                            #             best_errors[i] = error
                            #             best_gears[i] = (g1, g2, c1, c2, a, s)
                            #             best_ratios[i] = ratio
                            #             best_drive[i] = correction_drive
                            #     correction_drive = 's->p'
                            #     correction_ratio = 1.0 / ((s + a) / a)
                            #     ratio = g1 * g2 * (c2 / c1) * correction_ratio
                            #     for i,dr in enumerate(desired_ratios):
                            #         error = abs(abs(dr) - ratio)
                            #         if error < best_errors[i]:
                            #             best_errors[i] = error
                            #             best_gears[i] = (g1, g2, c1, c2, a, s)
                            #             best_ratios[i] = ratio
                            #             best_drive[i] = correction_drive


    return({'desired_ratios':desired_ratios, 'best_drive':best_drive, 'best_errors':best_errors, 'best_ratios':best_ratios, 'best_gears':best_gears})

def planetary_ok(g1, g2):
    # see if they're ok as planetary sun/ring
    if abs(g1 - g2) < 2 * 8:
        return False
    if (g1 ^ g2) & 1:
        return False
    return True

def master_gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic):
    """
    """
    gHr_range = range(8, 20+1)
    gMin_range = range(8, 20+1)
    g24a_range = range(16, 40+1) # This is the 24h rotator

    gMercG1_range   = range(1, 30+1) # Geneva stage 1
    gVenusG1_range  = range(1, 30+1) # Geneva stage 1
    gMarsG1_range   = range(1, 30+1) # Geneva stage 1
    gVenusG1_range  = range(1, 30+1) # Geneva stage 1
    gLunaG1_range   = range(1, 30+1) # Geneva stage 1

    gMerc1hub_range  = range(8, 40+1) # Connects to Hr or Min
    gVenus1hub_range = range(8, 40+1) # Connects to Hr or Min
    gMars1hub_range  = range(8, 40+1) # Connects to Hr or Min
    gVenus1hub_range = range(8, 40+1) # Connects to Hr or Min
    gLuna1hub_range  = range(8, 40+1) # Connects to Hr or Min

    gMerc1hubG_range  = range(8, 40+1) # Connects to Hr or Min
    gVenus1hubG_range = range(8, 40+1) # Connects to Hr or Min
    gMars1hubG_range  = range(8, 40+1) # Connects to Hr or Min
    gVenus1hubG_range = range(8, 40+1) # Connects to Hr or Min
    gLuna1hubG_range  = range(8, 40+1) # Connects to Hr or Min


    gMercG2_range   = range(1, 30+1) # Geneva stage 2
    gVenusG2_range  = range(1, 30+1) # Geneva stage 2
    gMarsG2_range   = range(1, 30+1) # Geneva stage 2
    gVenusG2_range  = range(1, 30+1) # Geneva stage 2
    gLunaG2_range   = range(1, 30+1) # Geneva stage 2

    gMercG2b_range   = range(8, 70+1) # teeth on G2
    gVenusG2b_range  = range(8, 70+1) # teeth on G2
    gMarsG2b_range   = range(8, 70+1) # teeth on G2
    gVenusG2b_range  = range(8, 70+1) # teeth on G2
    gLunaG2b_range   = range(8, 70+1) # teeth on G2

    gMercD_range   = range(8, 70+1) # teeth on driver
    gVenusD_range  = range(8, 70+1) # teeth on driver
    gMarsD_range   = range(8, 70+1) # teeth on driver
    gVenusD_range  = range(8, 70+1) # teeth on driver
    gLunaD_range   = range(8, 70+1) # teeth on driver

    gMercG1_range = [19]
    gMercG2_range = [20, 21, 22]
    debug_count = 0

    best_errors = [1.0] * len(desired_ratios)
    best_ratios = [1.0] * len(desired_ratios)
    best_gears = [None] * len(desired_ratios)
    best_drive = ['(none)'] * len(desired_ratios)
    do_planetaries = False
    iLuna = len(bodies) - 1
    iMerc = 0
    iVenus = 1
    iMars = 2
    for gMercG1 in gMercG1_range:
        for gMercG2 in gMercG2_range:
            for gMerc1hub in gMerc1hub_range:
                print_gear_ratio_search_progress(bodies, (-gMercG1, -gMercG2), best_errors, best_ratios, best_gears, is_sydonic)
                for gMerc1hubG in gMerc1hubG_range:
                    if planetary_ok(gMerc1hub, gMerc1hubG):
                        for gMercG2b in gMercG2b_range:
                            for gMercD in gMercD_range:
                                if planetary_ok(gMercG2b, gMercD):
                                    ratio = float(gMerc1hubG * gMercG1 * gMercG2 * gMercD) / float(gMerc1hub * gMercG2b)
                                    # print('{} / {} = {}'.format(gMerc1hubG * gMercG1 * gMercG2 * gMercD, gMerc1hub * gMercG2b, ratio))
                                    # debug_count += 1
                                    # if debug_count > 20:
                                    #     exit()
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (-gMercG1, gMerc1hub, gMerc1hubG, -gMercG2, gMercG2b, gMercD)
                                            best_ratios[i] = ratio
                                            best_drive[i] = ('s-a')
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (-gMercG1, gMerc1hub, gMerc1hubG, -gMercG2, gMercG2b, gMercD)
                                            best_ratios[i] = ratio
                                            best_drive[i] = ('s-a')
                                    for i,dr in enumerate(desired_ratios):
                                        error = abs(abs(dr) - ratio)
                                        if error < best_errors[i]:
                                            best_errors[i] = error
                                            best_gears[i] = (-gMercG1, gMerc1hub, gMerc1hubG, -gMercG2, gMercG2b, gMercD)
                                            best_ratios[i] = ratio
                                            best_drive[i] = ('s-a')


    return({'desired_ratios':desired_ratios, 'best_drive':best_drive, 'best_errors':best_errors, 'best_ratios':best_ratios, 'best_gears':best_gears})

def platter_gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic):
    """
    Modular platters

This is really great.
Early result with common shaft:

  searching (-19, 35, 19)
      best_errors [0.0013272135313400213, 0.0021241261645172926, 0.008511924298545637, 0.00038662585033222285]
      best_ratios [231.7530864197531, 1167.8433179723502, 1559.9035250463821, 59.06156462585033]
      best_gears [(-19, 35, 13, 27, 38, -6, 9, 35), (-19, 35, 13, 31, 38, -27, 8, 40), (-19, 35, 13, 11, 37, -23, 14, 40), (-19, 35, 13, 21, 37, -8, 32, 19)]
      error_degrees_per_year [1.5060627732894087, 0.47832736250932584, -1.435035582148089, 1.7215477380706343]      

Best for G1=19
desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
best_drive ['(none)', '(none)', '(none)', '(none)']
best_errors [0.0013272135313400213, 0.0021241261645172926, 0.008511924298545637, 0.00038662585033222285]
best_ratios [231.7530864197531, 1167.8433179723502, 1559.9035250463821, 59.06156462585033]
best_gears [(-19, 35, 13, 27, 38, -6, 9, 35), (-19, 35, 13, 31, 38, -27, 8, 40), (-19, 35, 13, 11, 37, -23, 14, 40), (-19, 35, 13, 21, 37, -8, 32, 19)]
bodies ['mercury', '  venus', '   mars', '   luna']
error_degrees_per_year [1.5060627732894087, 0.47832736250932584, -1.435035582148089, 1.7215477380706343]

Best for G1=20
desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
best_drive ['(none)', '(none)', '(none)', '(none)']
best_errors [0.0013272135313400213, 0.0021241261647446663, 0.006841079072728462, 5.424511094531681e-06]
best_ratios [231.7530864197531, 1167.84331797235, 1559.8881720430109, 59.06118342451109]
best_gears [(-20, 31, 13, 27, 19, -19, 15, 31), (-20, 31, 13, 10, 27, -19, 14, 38), (-20, 31, 13, 15, 37, -29, 10, 26), (-20, 31, 13, 27, 17, -17, 38, 25)]
bodies ['mercury', '  venus', '   mars', '   luna']
error_degrees_per_year [1.5060627732894087, 0.4783273625605276, -1.1533457706304298, 0.02415398452251853]

Best for G1=21
desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
best_drive ['(none)', '(none)', '(none)', '(none)']
best_errors [0.0008804843625966896, 0.009103356030664145, 0.003610407151199979, 9.141421059410959e-05]
best_ratios [231.75529411764703, 1167.8545454545454, 1559.8914027149324, 59.06126941421059]
best_gears [(-21, 34, 16, 25, 36, -6, 14, 38), (-21, 34, 16, 11, 31, -17, 15, 37), (-21, 34, 16, 13, 30, -18, 10, 38), (-21, 34, 16, 37, 28, -11, 39, 28)]
bodies ['mercury', '  venus', '   mars', '   luna']
error_degrees_per_year [0.9991344193360252, 2.0499649940145996, -0.6086828954645581, 0.4070445039838193]

Best for G1=22
desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
best_drive ['(none)', '(none)', '(none)', '(none)']
best_errors [0.0013272135313115996, 0.0013875176180135895, 0.003346455416703975, 0.00021227236759813195]
best_ratios [231.75308641975312, 1167.8440545808967, 1559.8916666666669, 59.0609657276324]
best_gears [(-22, 27, 11, 33, 38, -13, 22, 38), (-22, 27, 11, 11, 32, -23, 19, 37), (-22, 27, 11, 10, 39, -21, 16, 34), (-22, 27, 11, 21, 16, -8, 37, 40)]
bodies ['mercury', '  venus', '   mars', '   luna']
error_degrees_per_year [1.506062773257157, 0.31245208205911135, -0.564182954242549, 0.9451955009719262]

Best for G1=23:
desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
best_drive ['(none)', '(none)', '(none)', '(none)']
best_errors [0.001196522325699334, 0.0011141351380956621, 0.0014832285537522694, 0.0001715988620176745]
best_ratios [231.75561015561013, 1167.8465562336528, 1559.8935298935298, 59.06100640113798]
best_gears [(-23, 37, 13, 33, 26, -13, 10, 28), (-23, 37, 13, 31, 32, -28, 8, 40), (-23, 37, 13, 11, 35, -26, 15, 35), (-23, 37, 13, 38, 21, -23, 40, 23)]
bodies ['mercury', '  venus', '   mars', '   luna']
error_degrees_per_year [1.3577602168703062, 0.25088967453370076, -0.2500592905244979, 0.7640866033871646]

24-30:
      best_errors [0.0002094545855868546, 0.002124126164289919, 0.0020335698811777547, 0.00014134032982582312]
      best_ratios [231.75420417869884, 1167.8433179723504, 1559.8970466919648, 59.061319340329824]
      best_gears [(-27, 29, 13, 29, 23, -13, 21, 39), (-27, 29, 13, 31, 38, -29, 14, 38), (-27, 29, 13, 11, 26, -28, 19, 37), (-27, 29, 13, 23, 13, -7, 30, 37)]
      error_degrees_per_year [0.2376797302003289, 0.4783273624581241, -0.34284199857996345, 0.6293529646316143]

    """
    gG1_range = range(6, 30+1)
    gG1r_range = range(8, 40+1)
    gShaft_range = range(8, 40+1)

    gP1a_range = range(8, 40+1)
    gP1b_range = range(8, 40+1)


    gGP_range = range(6, 30+1)
    gPDa_range = range(8, 40+1)
    gPDb_range = range(8, 40+1)


    gG1_range = [18]
    # gG1r_range = [35]
    # gShaft_range = [17]

    debug_count = 0

    current_errors = [1.0] * len(desired_ratios)
    current_gears = [1.0] * len(desired_ratios)
    best_max_error_normalized = 1000000.0
    best_errors = [1.0] * len(desired_ratios)
    best_ratios = [1.0] * len(desired_ratios)
    best_gears = [None] * len(desired_ratios)
    best_drive = ['(none)'] * len(desired_ratios)
    do_planetaries = False
    iLuna = len(bodies) - 1
    iMerc = 0
    iVenus = 1
    iMars = 2
    require_planetary = (1, 1, 1, 0)
    for gG1 in gG1_range:
        for gG1r in gG1r_range:
            for gShaft in gShaft_range:
                if planetary_ok(gG1r, gShaft) and gG1r > gShaft:
                    print_gear_ratio_search_progress(bodies, (-gG1, gG1r, gShaft), best_errors, best_ratios, best_gears, is_sydonic)
                    shaft_ratio = float(gG1 * gShaft) / float(gG1r)
                    bodies_ok = True
                    max_body_error_normalized = 0.0
                    new_best_gears = [None] * len(desired_ratios)
                    new_best_errors = [None] * len(desired_ratios)
                    new_best_ratios = [None] * len(desired_ratios)
                    # Check all the bodies. If the max normalized ratio is worse, don't take it.
                    for body_index,desired_ratio in enumerate(desired_ratios):
                        desired_ratio = abs(desired_ratio)
                        best_body_error = 1000000.0
                        best_body_error_normalized = 1000000.0
                        for gGP in gGP_range:
                            for gP1a in gP1a_range:
                                for gP1b in gP1b_range:
                                    for gPDa in gPDa_range:
                                        for gPDb in gPDb_range:
                                            if (planetary_ok(gPDa, gPDb) and gPDb > gPDa) or not require_planetary[body_index]:
                                                ratio = float(shaft_ratio * gP1b * gGP * gPDb) / float(gP1a * gPDa)
                                                error = abs(desired_ratio - ratio)
                                                if error < best_body_error:
                                                    best_body_gears = (-gG1, gG1r, gShaft, gP1a, gP1b, -gGP, gPDa, gPDb)
                                                    best_body_error_normalized = error / desired_ratio
                                                    best_body_error = error
                                                    best_body_ratio = ratio
                        if best_body_error_normalized > max_body_error_normalized:
                            max_body_error_normalized = best_body_error_normalized
                        if max_body_error_normalized >= best_max_error_normalized:
                            bodies_ok = False
                            break
                        new_best_gears[body_index] = best_body_gears
                        new_best_errors[body_index] = best_body_error
                        new_best_ratios[body_index] = best_body_ratio
                    if bodies_ok:
                        best_max_error_normalized = max_body_error_normalized
                        best_gears = new_best_gears
                        best_errors = new_best_errors
                        best_ratios = new_best_ratios


    return({'desired_ratios':desired_ratios, 'best_drive':best_drive, 'best_errors':best_errors, 'best_ratios':best_ratios, 'best_gears':best_gears})


def do_abs_search():
# desired_ratios [175.9382, 449.402, 730.512726008, 1373.942]
# best_drive ['p->s', 'a->s', 'a->p', 'a->p']
# best_errors [4.561403505931594e-05, 8.612440183242143e-05, 9.450482036754693e-05, 0.00020105820090066118]
# best_ratios [175.93824561403505, 449.40191387559815, 730.5128205128203, 1373.941798941799]
# best_gears [(22, 11, 19, 28, 38, 37), (17, 13, 11, 17, 25, 19), (15, 11, 39, 37, 33, 9), (17, 13, 9, 25, 26, 21)]
# bodies ['mercury', '  venus', '  earth', '   mars']
# error_degrees_per_year [0.06818182699576289, 0.05039895185385869, 0.034021735332316894, 0.03848430777531672]
    chain_name = 'chain-v1'
    bodies = planets
    desired_ratios = [b['em_gear_ratio'] for b in bodies]
    result = gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic=False)
    result['bodies'] = [b['nspace']+b['name'] for b in bodies]
    result['error_degrees_per_year'] = [None]*len(bodies)

    for i,body in enumerate(bodies):
        error = result['best_errors'][i]
        days_per_orbit = body['orbital_period_days']
        earth_years_per_orbit = days_per_orbit / earth['orbital_period_days']
        orbit_per_earth_year = 1.0 / earth_years_per_orbit
        result['error_degrees_per_year'][i] = error * orbit_per_earth_year * 360.0

    for key,val in result.items():
        print(key, val)
    return result

def do_sydonic_search():
# desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
# best_drive ['s->p', 'p->a', 'a->p', 's->p']
# best_errors [5.269745307145968e-05, 9.361577099298302e-05, 0.00027627997815216077, 0.00011922024233257389]
# best_ratios [231.7544663307375, 1167.8455357142857, 1559.8947368421054, 59.06129722024233]
# best_gears [(24, 17, 37, 31, 40, 19), (29, 23, 28, 37, 40, 13), (22, 21, 19, 29, 40, 33), (9, 9, 23, 33, 31, 30)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [0.05979872148979634, 0.021081132362270317, -0.046578374686813084, 0.5308577746250345]

# Using 22 as the drive Geneva is pretty good (and then 11 does both Mercury ans Venus):
# best_drive ['p->s', 'a->p', 'a->p', 's->p']
# best_errors [0.0003813752198880138, 0.000280808191973847, 0.00027627997815216077, 0.00014513615023048487]
# best_ratios [231.75403225806454, 1167.8451612903227, 1559.8947368421054, 59.06103286384977]
# best_gears [(22, 11, 8, 19, 37, 25), (22, 11, 10, 22, 37, 31), (22, 21, 19, 29, 40, 33), (22, 17, 33, 10, 37, 34)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [0.43276760503525186, 0.06323458751254847, -0.046578374686813084, 0.6462547988627]

# This is with NO planetary gears!
# This is pretty good for no planetaries.
# Time to search each planet separately, and Luna clearly needs her own base setup
# desired_ratios [231.75441363328443, 1167.8454420985147, -1559.8950131220836, 59.061178]
# best_drive ['a->s', 'a->s', 'a->s', 'a->s']
# best_errors [0.0003813752198880138, 0.000280808191973847, 0.003610407151199979, 0.0005503950617224973]
# best_ratios [231.75403225806454, 1167.8451612903227, 1559.8914027149324, 59.06172839506172]
# best_gears [(19, 11, 8, 25, 11, 31), (22, 22, 10, 22, 34, 31), (21, 16, 13, 27, 38, 17), (4, 2, 9, 23, 26, 9)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [0.43276760503525186, 0.06323458751254847, -0.6086828954645581, 2.450770875096454]

# At 19, no planetaries:
# This is the best so far, but Luna's gears are huge, she probably needs her own train.
# best_gears [(19, 11, 8, 25, 11, 31), (19, 20, 8, 31, 23, 29), (19, 21, 17, 24, 36, 13), (19, 1, 8, 27, 35, 38)]
# error_degrees_per_year [0.43276760503525186, 0.13838069333711012, -0.6086828954645581, 5.8865337322294256]

# at 22 no planetaries:
# best_gears [(22, 19, 16, 25, 11, 31), (22, 22, 10, 22, 34, 31), (22, 25, 14, 27, 25, 17), (22, 23, 33, 13, 8, 27)]
# bodies ['mercury', '  venus', '   mars', '   luna']
# error_degrees_per_year [0.43276760503525186, 0.06323458751254847, -3.5511279889928047, 2.450770875096454]

    chain_name = 'chain-v2'
    bodies = planets_without_earth + moons
    desired_ratios = [b['es_gear_ratio'] for b in bodies]
    result = platter_gear_ratio_search(chain_name, bodies, desired_ratios, is_sydonic=True)
    result['bodies'] = [b['nspace']+b['name'] for b in bodies]
    result['error_degrees_per_year'] = [None]*len(bodies)

    for i,body in enumerate(bodies):
        error = result['best_errors'][i]
        days_per_relative_orbit = body['sydonic_period_days']
        earth_years_per_relative_orbit = days_per_relative_orbit / earth['orbital_period_days']
        relative_orbit_per_earth_year = 1.0 / earth_years_per_relative_orbit
        result['error_degrees_per_year'][i] = error * relative_orbit_per_earth_year * 360.0

    for key,val in result.items():
        print(key, val)
    return result

"""
Favorite results:

# bodies ['mercury', '  venus', '   mars', '   luna']
# negative teeth counts indicate Geneva gears

# Thoughts:
1. G should not be last, as step size will be too big. Last should be a step-down.
2. Maybe the G->G->c->c can be part-reduced by putting one c as the G link, and keep it concentric

G->G->corrrection->correction (total planet parts with shared G1: 1+1+3+3+3=11)
Everything but Luna:
    # best_gears [(-19, -11, 8, 25, 11, 31), (-19, -20, 8, 31, 23, 29), (-19, -21, 17, 24, 36, 13), (-19, -1, 8, 27, 35, 38)]
    # error_degrees_per_year [0.43276760503525186, 0.13838069333711012, -0.6086828954645581, 5.8865337322294256]
Luna options:
    # best_gears [(-19, -11, 8, 25, 11, 31), (-22, -22, 10, 22, 34, 31), (-21, -16, 13, 27, 38, 17), (-4, -2, 9, 23, 26, 9)]
    # error_degrees_per_year [0.43276760503525186, 0.06323458751254847, -0.6086828954645581, 2.450770875096454]
    # best_gears [(-22, -19, 16, 25, 11, 31), (-22, -22, 10, 22, 34, 31), (-22, -25, 14, 27, 25, 17), (-22, -23, 33, 13, 8, 27)]
    # error_degrees_per_year [0.43276760503525186, 0.06323458751254847, -3.5511279889928047, 2.450770875096454]

Snowman->G->G->correction: (total planet parts with shared snowman core: 1+3+3+3+3=13)
All driven from 9:
    best_gears [(9, 8, -13, -19, 18, 19), (9, 11, -7, -28, 8, 39), (9, 31, -12, -28, 23, 31), (9, 8, -1, -23, 9, 26)]
    bodies ['mercury', '  venus', '   mars', '   luna']
    error_degrees_per_year [1.5060627732894087, 2.7267465601781287, -1.8469421297166195, 2.4507708751280926]

1:2->correction->G->G->corrrection (total planet parts with shared snowman core: 1+1+3+3+3+3=14)
    # best_gears [(10, 20, 31, 11, -19, -25, 11, 16), (10, 20, 10, 22, -17, -22, 22, 31), (10, 20, 17, 37, -21, -23, 23, 31), (10, 20, 27, 13, -23, 9, 24)]
    # error_degrees_per_year [0.43276760503525186, 0.06323458746134669, -0.42486412657231726, 2.450770875064815]

"""

#def eval_master():



if __name__ == "__main__":
    list_planets()
    make_conjunctions()
    make_gear_ratios()
#    do_abs_search()
    do_sydonic_search()
