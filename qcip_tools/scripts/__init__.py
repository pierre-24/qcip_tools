"""
Scripts to ease the manipulation of quantum chemistry results in Python 3
"""

import sys


def exit_failure(msg, status=1):
    """Write a message in stderr and exits

    :param msg: the msg
    :type msg: str
    :param status: exit status (!=0)
    :type status: int
    """

    sys.stderr.write(msg)
    sys.stderr.write('\n')
    return sys.exit(status)
