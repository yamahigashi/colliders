"""Run with mayapy; one plugin and one standalone session per invocation."""

from pathlib import Path
import sys
import unittest

from helpers import maya_session, parser


def main():
    args_parser = parser(__doc__)
    args_parser.add_argument("--pattern", default="test_*.py")
    args = args_parser.parse_args()
    with maya_session(args.plugin):
        suite = unittest.defaultTestLoader.discover(str(Path(__file__).parent), args.pattern)
        if suite.countTestCases() == 0:
            raise RuntimeError("No tests discovered")
        result = unittest.TextTestRunner(verbosity=2).run(suite)
        return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    sys.exit(main())
