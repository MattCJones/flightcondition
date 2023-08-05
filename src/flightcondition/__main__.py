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
SHORT_DESCRIPTION = ("Compute airspeed (true/calibrated/equivalent/Mach), "
                     "atmospheric data, and other flight condition quantities,"
                     " with easy unit conversion.")


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
        '--h', dest='altitude', metavar='', nargs=2, type=str,
        default=[0, 'kft'],
        help=("altitude and length unit, default='0 kft'; aliases include: "
              "--altitude --Altitude"))
    parser.add_argument(
        '--altitude', dest='altitude', metavar='', nargs=2,
        type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--Altitude', dest='altitude', metavar='', nargs=2,
        type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--kft', dest='kft', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--km', dest='km', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--M', dest='M', metavar='', nargs=None, type=float, default=None,
        help=("Mach number, e.g. '0.5'; aliases include: "
              "--mach-number --Mach"))
    parser.add_argument(
        '--mach-number', dest='M', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--Mach', dest='M', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--TAS', dest='TAS', metavar='', nargs=2, type=str, default=None,
        help="true airspeed and speed unit, e.g. '150 knots'")
    parser.add_argument(
        '--tas', dest='TAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--KTAS', dest='KTAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--ktas', dest='KTAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--CAS', dest='CAS', metavar='', nargs=2, type=str, default=None,
        help="calibrated airspeed and speed unit, e.g. '150 knots'")
    parser.add_argument(
        '--cas', dest='CAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--KCAS', dest='KCAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--kcas', dest='KCAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--EAS', dest='EAS', metavar='', nargs=2, type=str, default=None,
        help="equivalent airspeed and speed unit, e.g. '150 knots'")
    parser.add_argument(
        '--eas', dest='EAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--KEAS', dest='KEAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--keas', dest='KEAS', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--L', dest='length_scale', metavar='', nargs=2,
        type=str, default=None,
        help=("length scale, e.g. '10 ft'; aliases include: "
              "--length-scale --Length --ell"))
    parser.add_argument(
        '--length-scale', dest='length_scale', metavar='', nargs=2,
        type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--Length', dest='length_scale', metavar='', nargs=2,
        type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--ell', dest='length_scale', metavar='', nargs=2,
        type=str, default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--ft', dest='ft', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--m', dest='m', metavar='', nargs=None, type=float,
        default=None, help=argparse.SUPPRESS)
    parser.add_argument(
        '--units', dest='units', metavar='', nargs=None,
        type=str, default="", help=f"Unit system, i.e. {dir(unit.sys)}")
    parser.add_argument(
        '--no-full-output', dest='no_full_output', default=False,
        action='store_true', help="display abbreviated output")
    parser.add_argument(  # hidden option to turn off pretty print
        '--no-pretty-print', dest='no_pretty_print', default=False,
        action='store_true', help=argparse.SUPPRESS)

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
    h_ = h_ if args.kft is None else args.kft*unit('kft')
    h_ = h_ if args.km is None else args.km*unit('km')
    M_ = None if args.M is None else args.M
    TAS_ = None if args.TAS is None else _dimension(args.TAS)
    CAS_ = None if args.CAS is None else _dimension(args.CAS)
    EAS_ = None if args.EAS is None else _dimension(args.EAS)
    TAS_ = TAS_ if args.KTAS is None else args.KTAS*unit('knots')
    CAS_ = CAS_ if args.KCAS is None else args.KCAS*unit('knots')
    EAS_ = EAS_ if args.KEAS is None else args.KEAS*unit('knots')
    L_ = None if args.length_scale is None else _dimension(args.length_scale)
    L_ = L_ if args.ft is None else args.ft*unit('ft')
    L_ = L_ if args.m is None else args.m*unit('m')
    fc = FlightCondition(h=h_, M=M_, TAS=TAS_, CAS=CAS_, EAS=EAS_, L=L_,
                         units=args.units)

    # Print output but catch common unicode exception
    full_output = not args.no_full_output
    try:
        print(fc.tostring(full_output=full_output,
              pretty_print=(not args.no_pretty_print)))
    except UnicodeEncodeError:
        print(fc.tostring(full_output=full_output, pretty_print=False))


if __name__ == '__main__':
    main()
