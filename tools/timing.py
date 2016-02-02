#!/usr/bin/env python

"""
Usage: ./timing.py <number_of_measures> <executable1> <executable2> <executable3> ... 
Executables have to print time of their execution in microseconds 
at the end of their output after whitespace 
"""


def run(args):
    import subprocess
    popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    popen.wait()
    return popen.stdout.read()


def timing(times, exes):
    for trying in range(len(times[0])):
        for exe_index in range(len(times)):
            print "Run " + exes[exe_index],
            time = int(run([exes[exe_index]]).split(" ")[-1])
            print "Time = " + str(time)
            times[exe_index][trying] = float(time)
    return times


def square_product(la, lb):
    plain_product = [float(a * b) for a, b in zip(la, lb)]
    return sum(plain_product)


def average(times):
    from math import sqrt
    mean = [sum(e) / len(e) for e in times]
    sq_mean = [square_product(e, e) / len(e) for e in times]
    stdev = [sqrt(sq_mean_ - mean_*mean_) for sq_mean_, mean_ in zip(sq_mean, mean)]
    return [[m, e] for m, e in zip(mean, stdev)]


def main():
    print
    import sys
    args = sys.argv
    number_of_tryings = int(sys.argv[1])
    exes = args[2:]
    number_of_exes = len(exes)
    assert number_of_exes > 0, 'specify executables'

    times = [[0.0 for i in range(number_of_tryings)] for i in range(number_of_exes)]
    times = timing(times, exes)

    av = average(times)

    print
    for i in range(number_of_exes):
        print exes[i] + " : " + str(av[i][0]) + " +/- " + str(av[i][1])

    print
    base = min([a[0] for a in av])
    for i in range(number_of_exes):
        print exes[i] + " : " + str(av[i][0] / base) + " +/- " + str(av[i][1] / base)


main()
