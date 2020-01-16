#!/usr/bin/env python
# -*- coding: utf-8 -*-
''' Generate IN3D.out and IN3D_gmrot.out
    Usage: python generate_in3dout.py 
           param.sh should contain the assignment of variables mentioned below
'''


import os
import argparse

class CustomArgumentParser(argparse.ArgumentParser):
    def __init__(self, *args, **kwargs):
        super(CustomArgumentParser, self).__init__(*args, **kwargs)

    def convert_arg_line_to_args(self, line):
        for arg in line.split():
            if not arg.strip():
                continue
            if arg[0] == '#':
                break
            yield arg

def convert_arg_line_to_args(arg_line):
    for arg in arg_line.split():
        if not arg.strip():
            continue
        if arg[0] == '#':
            break
        yield arg

def generate_files():
    # Replace the method of the argparse.ArgumentParser class
    # So that we don't have to build a CustomArgumentParser class.
    parser = argparse.ArgumentParser(description="Check if param.sh and IN3D.out (for PGV2.lsf and extrts.lsf) and IN3D_gmrot.out (for gmrot.lsf) are consistent", fromfile_prefix_chars="@")
    parser.convert_arg_line_to_args = convert_arg_line_to_args

    ### generate a child class inheritating the parant argparse.ArgumentParser class, with a modified convert_arg_line_to_args method. Equivalently working
    # parser = CustomArgumentParser(argparse.ArgumentParser(), description="Check if param.sh and IN3D.out (for PGV2.lsf and extrts.lsf) and IN3D_gmrot.out (for gmrot.lsf) are consistent", fromfile_prefix_chars="@")

    parser.add_argument('--NBGX', type=int)
    parser.add_argument('--NBGY', type=int)
    parser.add_argument('--NBGZ', type=int)
    parser.add_argument('--NEDX', type=int)
    parser.add_argument('--NEDY', type=int)
    parser.add_argument('--NEDZ', type=int)
    parser.add_argument('--NSKPX', type=int, default=1)
    parser.add_argument('--NSKPY', type=int, default=1)
    parser.add_argument('--NSKPZ', type=int, default=1)
    parser.add_argument('--DT', type=float)
    parser.add_argument('--TMAX', type=float)
    parser.add_argument('--NTISKP', type=int, default=10)
    parser.add_argument('--SXRGO', default=repr('output_sfc/SX_0_'))
    parser.add_argument('--SYRGO', default=repr('output_sfc/SY_0_'))
    parser.add_argument('--SZRGO', default=repr('output_sfc/SZ_0_'))
    parser.add_argument('--READ_STEP', type=int, default=100)
    parser.add_argument('--WRITE_STEP', type=int, default=100)
    parser.add_argument('--IVELOCITY', type=int, default=1)

    args = parser.parse_known_args(['@param.sh'])[0]
    for key, value in vars(args).items():
        print(key, value, '\n')

    with open('IN3D_gmrot.out', 'w') as f_in3dgmrot, open('IN3D.out', 'w') as f_in3dout:
        for key, value in vars(args).items():
            f_in3dout.write(f'{value} {key}\n')
            if key == 'WRITE_STEP' and vars(args)['IVELOCITY'] == 0:
                write_gmrot = int(int(args.TMAX / args.DT) / args.NTISKP)
                f_in3dgmrot.write(f'{write_gmrot} {key}\n')
                continue
            f_in3dgmrot.write(f'{value} {key}\n')


if __name__ == "__main__":
    generate_files()
