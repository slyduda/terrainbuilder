from simplexnoise import simplexnoise

SQRDMD_ROUGHNESS = 0.35


def nearest_pow(num):
    n = 1
    while n < num:
        n << 1
    return n


def generate_noise(tmp, tmp_dim, randstate, use_simplex):
    if use_simplex:
        simplexnoise(randstate.next(), tmp, tmp_dim.get_width(),
                     tmp_dim.get_height(), SQRDMD_ROUGHNESS)
    else:
        side = tmp_dim.get_max()
        side = nearest_pow(side)+1
        square_tmp = [0.0 for _ in range(side*side)]
        for y in range(tmp_dim.get_height()):
            for x in range(tmp_dim.get_width(), side):
                # we simply put it as a mix between the east and west border(they should be fairly
                # similar because it is a toroidal world)
                square_tmp[y*side+x] = (square_tmp[y*side+0] +
                                        square_tmp[y*side+(tmp_dim.get_width()-1)])/2

        # 2) below the valid area
        for y in range(tmp_dim.get_height(), side):
            for x in range(side):
                # we simply put it as a mix between the north and south border(they should be fairly
                # similar because it is a toroidal world)
                square_tmp[y*side+x] = (square_tmp[(0)*side+x] +
                                        square_tmp[(tmp_dim.get_height()-1)*side+x])/2

        sqrdmd(randstate.next(), square_tmp, side, SQRDMD_ROUGHNESS)

        # Calcuate deltas(noise introduced)
        deltas = [0.0 for _ in tmp_dim.get_width()*tmp_dim.get_height()]
        for y in range(tmp_dim.get_height()):
            for x in range(tmp_dim.get_width()):
                deltas[y*tmp_dim.get_width()+x] = square_tmp[y*side+x] - \
                    tmp[y*tmp_dim.get_width()+x]

        # make it tileable
        for y in range(tmp_dim.get_height()):
            for x in range(tmp_dim.get_width()):
                specular_x = tmp_dim.get_width() - 1 - x
                specular_y = tmp_dim.get_height() - 1 - y
                my_delta = deltas[y * tmp_dim.get_width() + x]
                specular_width_delta = deltas[y *
                                              tmp_dim.get_width() + specular_x]
                specular_height_delta = deltas[specular_y *
                                               tmp_dim.get_width() + x]
                opposite_delta = deltas[specular_y *
                                        tmp_dim.get_width() + specular_x]
                tmp[y * tmp_dim.get_width() + x] += (my_delta + specular_width_delta +
                                                     specular_height_delta + opposite_delta) / 4
