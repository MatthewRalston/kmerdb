#!/usr/bin/env python


logger = None


if __name__ == '__main__':
    from kdb import cli

    if len(sys.argv) == 1:
        sys.exit(-1)
    else:
        cli()
