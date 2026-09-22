import sys, os, time
from signal import SIGTERM

"""
    So I wrote this daemon to keep track of long=running simulations, which run
    on SLURM, which I want to inherit.

    The recipe for building this daemon is from

    References:
        Fork a Daemon Process on UNIX,
            Hermann, Jurgen.  2001.
            http://code.activestate.com/recipes/66012/#c9

        UNIX Programming FAQ
            1.7 How do I get my program to act like a daemon?
                http://www.erlenstar.demon.co.uk/unix/faq_2.html#SEC16

        Advanced Programming in the Unix Environment
            W. Richard Stevens, 1992, Addison-Wesley, ISBN 0-201-56317-7.

    History:

"""
def main(path):

    f = open(path, 'w')
    while True:
        f.write('{}\n'.format(time.ctime(time.time())))
        f.flush()
        time.sleep(sleep_time)

if __name__ == "__main__":

    # UNIX double fork magic
    # do the UNIX double-fork magic, see Stevens' "Advanced
    # Programming in the UNIX Environment" for details (ISBN 0201563177)
    try:
        pid = os.fork()
        if pid > 0:
            sys.exit(0)
    except OSError, e:
        print >>sys.stderr, "fork #1 failed: %d (%s)" % (e.errno, e.strerror)
        sys.exit(0)

    os.chdir("/")
    os.setid()
    os.umask(0)

    try:
        pid = os.fork()
        if pid > 0:
            # exit from second parent, print eventual PID score before
