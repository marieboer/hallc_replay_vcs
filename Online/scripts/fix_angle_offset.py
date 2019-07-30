#!/usr/bin/env python3

from epics import caget, caput
import argparse

class InvalidArgumentError(Exception):
    pass

parser = argparse.ArgumentParser()
parser.add_argument(
        "-s", "--spectrometer",
        required=True,
        help="Spectrometer arm (HMS or SHMS)",
        dest="spectrometer")
parser.add_argument(
        "-a", "--angle",
        required=True,
        help="Correct spectrometer angle",
        dest="angle")
args = parser.parse_args()

if args.spectrometer not in ['HMS', 'SHMS']:
    print("Unknown value for spectrometer", args.spectrometer)
    args.print_help()
    raise InvalidArgumentError()

true_angle = float(args.angle)
print("True {} angle: ".format(args.spectrometer), true_angle)
if true_angle <= 3:
    print("Invalid spectrometer angle", args.spectrometer)
    raise InvalidArgumentError()

epics_angle = caget("ec{}_Angle".format(args.spectrometer))
print("Uncorrected angle from EPICS:", epics_angle)

offset = round(true_angle - epics_angle, 4)
print("Calculated angle offset: ", offset)
caput("hc{}AngleEncoderOffset".format(args.spectrometer), offset)
print("EPICS corrected angle: ", epics_angle +
        caget("hc{}AngleEncoderOffset".format(args.spectrometer)))
