#Modified by sam ward the fast ion legend
#inherited from ben dudson and nick walkden

def dump_equilibrium_GEQDSK(output_data, filepath):
    """
    generic function for writing G-EQDSK-formatted data to file

    notes:
        originally written by Ben Dudson and edited by Nick Walkden
    """

    def write_number(file, number, counter):

        if number == 0:
            separator = "  "
            number = np.abs(number)
        elif number < 0:
            separator = " -"
            number = np.abs(number)
        else:
            separator = "  "
        if utils.get_next(counter) == 4:
            last = "\n"
        else:
            last = ""

        string = '{:.8e}'.format(float(number))
        # mant,exp = string.split('E')
        file.write(separator + string + last)

    def write_1d(file, array, counter):
        for num in array:
            write_number(file, num, counter)

    def write_2d(file, array, counter):
        ny = array.shape[1]
        for j in np.arange(ny):
            write_1d(file, array[:, j], counter)

    def write_bndry(file, R, Z, counter):
        for i in np.arange(len(list(R))):
            write_number(file, R[i], counter)
            write_number(file, Z[i], counter)

    print("writing equilibrium to GEQDSK")

    with open(filepath, 'w') as file:

        cnt = itertools.cycle([0, 1, 2, 3, 4])  # initialise counter

        EFIT_shot = 19113  # just 'make up' a shot number and time (in ms) for now
        EFIT_time = 23
        line = "LOCUSTIO   " + time.strftime("%d/%m/%y") + "      #" + utils.fortran_string(EFIT_shot,
                                                                                            6) + utils.fortran_string(
            EFIT_time, 6) + utils.fortran_string(output_data['xdum'], 14) + utils.fortran_string(output_data['nR_1D'],
                                                                                                 4) + utils.fortran_string(
            output_data['nZ_1D'], 4) + "\n"
        file.write(line)

        float_keys = [
            'rdim', 'zdim', 'rcentr', 'rleft', 'zmid',
            'rmaxis', 'zmaxis', 'simag', 'sibry', 'bcentr',
            'current', 'simag', 'xdum', 'rmaxis', 'xdum',
            'zmaxis', 'xdum', 'sibry', 'xdum', 'xdum']
        for key in float_keys:
            write_number(file, output_data[key], cnt)

        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset the counter (otherwise our newlines will be out of sync)
        write_1d(file, output_data['fpol'], cnt)

        file.write("\n")
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_1d(file, output_data['pres'], cnt)

        file.write("\n")
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_1d(file, output_data['ffprime'], cnt)

        file.write("\n")
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_1d(file, output_data['pprime'], cnt)

        file.write("\n")
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_2d(file, output_data['psirz'], cnt)

        file.write("\n")
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_1d(file, output_data['qpsi'], cnt)

        file.write(
            "\n" + utils.fortran_string(len(output_data['lcfs_r']), 5) + utils.fortran_string(len(output_data['rlim']),
                                                                                              5))  # write out number of limiter/plasma boundary points

        file.write("\n")
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_bndry(file, output_data['lcfs_r'], output_data['lcfs_z'], cnt)

        file.write("\n")  # file has a newline here
        cnt = itertools.cycle([0, 1, 2, 3, 4])  # reset again
        write_bndry(file, output_data['rlim'], output_data['zlim'], cnt)


print("finished writing equilibrium to GEQDSK")
