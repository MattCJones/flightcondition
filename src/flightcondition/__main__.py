"""Command Line Interface (CLI) main function for easy use from the terminal.

Dependencies: numpy, pint

Author: Matthew C. Jones
Email: matt.c.jones.aoe@gmail.com

:copyright: 2021 Matthew C. Jones
:license: MIT License, see LICENSE for more details.
"""

import argparse
import sys

from .flightcondition import FlightCondition
from .units import unit

RUNCMD = 'flightcondition'
SHORT_DESCRIPTION = ("Airspeed conversions (true/calibrated/equivalent/Mach), "
                     "atmospheric data, and more with built-in unit checking.")


def _parse_args(args_arr):
    """Parse arguments.

    Args:
        args (Sequence[str]): argument array

    Returns:
        argparse.Namespace: parsed arguments array
    """
    # Initialize parser and define arguments
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter, prog=RUNCMD,
        description=SHORT_DESCRIPTION)
    parser.add_argument(
        '--altitude', dest='altitude', metavar='', nargs=2, type=str,
        default=[0, 'kft'],
        help="altitude and lenght unit, default='0 kft'")
    parser.add_argument(
        '-M', '--mach-number', dest='M', metavar='', nargs=None, type=float,
        default=None, help="Mach number")
    parser.add_argument(
        '-T', '--TAS', dest='TAS', metavar='', nargs=2, type=str, default=None,
        help="true airspeed and speed unit, e.g. '150 knots'")
    parser.add_argument(
        '-C', '--CAS', dest='CAS', metavar='', nargs=2, type=str, default=None,
        help="calibrated airspeed and speed unit, e.g. '150 knots'")
    parser.add_argument(
        '-E', '--EAS', dest='EAS', metavar='', nargs=2, type=str, default=None,
        help="equivalent airspeed and speed unit, e.g. '150 knots'")
    parser.add_argument(
        '-L', '--length-scale', dest='length_scale', metavar='', nargs=2,
        type=str, default=None,
        help="length scale, e.g. '10 ft'")
    parser.add_argument(
        '--no-pretty-print', dest='no_pretty_print', default=False,
        action='store_true',
        help="turn off pretty print (throws error in some terminals)")

    args = parser.parse_args(args_arr)

    return args


def _dimension(num_unit_arr):
    """Parse length 2 array into dimensional unit quantity.

    Args:
        array: float value and unit string

    Returns:
        dimensional: dimensionalized quantity
    """
    val, unit_ = num_unit_arr
    return float(val) * unit(unit_)


def main():
    """Main function. """
    args_arr = sys.argv[1:] if len(sys.argv) > 1 else None
    args = _parse_args(args_arr)
    h_ = _dimension(args.altitude)
    M_ = None if args.M is None else args.M
    TAS_ = None if args.TAS is None else _dimension(args.TAS)
    CAS_ = None if args.CAS is None else _dimension(args.CAS)
    EAS_ = None if args.EAS is None else _dimension(args.EAS)
    L_ = None if args.length_scale is None else _dimension(args.length_scale)
    fc = FlightCondition(h=h_, M=M_, TAS=TAS_, CAS=CAS_, EAS=EAS_, L=L_)
    print(fc.tostring(full_output=True,
                      pretty_print=(not args.no_pretty_print)))


if __name__ == '__main__':
    main()
