def sortPatientbyCancer(patientsCancer):
    with open(patientsCancer) as file:
        line = file.readline()
        line2 = file.readline()

    line2 = line2.split()[1:]

    from collections import Counter
    Counter(line2).keys()

    Types = Counter(line2).keys()
    TypesDiction = dict((k, []) for k in Types)

    for kind in TypesDiction:
        cancerind = [i for i, x in enumerate(line2) if x == kind]
        TypesDiction[kind] = cancerind

    return TypesDiction

