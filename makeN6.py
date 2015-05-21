class OptionError(ValueError):
    pass


def convert_repr_int(input_int, target_repr, digit):

    '''
    Return the target_repr representation of an integer
    '''

    if input_int > (target_repr ** digit - 1):
        print 'Cannot change representation'
        raise OptionError()

    output = str('')
    k = input_int

    for i in reversed(range(digit)):
        j, k = divmod(k, target_repr ** i)
        output += str(j)

    return output


def replaces_str(input_str, dict):
    output = input_str
    for key in dict.keys():
        output = output.replace(key, dict[key])
    return output


def make_random_oligo(N):
    oligo_quad = [convert_repr_int(i, 4, 6) for i in range(4 ** N)]
    dict_gatc = {'0': 'G', '1': 'A', '2': 'T', '3': 'C'}
    oligo_gatc = [replaces_str(s, dict_gatc) for s in oligo_quad]
    return oligo_gatc


if __name__ == '__main__':

    test = make_random_oligo(6)
