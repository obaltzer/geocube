#!/usr/bin/env python

import getopt
import sys
import os
import random
import signal

filename = None

def int_handler(signum, frame):
    print "Caught SIGINT. Deleting output file."
    if filename != None:
        os.remove(filename)
    sys.exit(123)
        
def usage():
    print """
    datagen.py [--dimf n] [--dimz=d1,d2,d3...] [--prec n] 
               [--number num] [filename]
    """
    
signal.signal(signal.SIGINT, int_handler)

short_opts = ""
long_opts = ['dimf=', 'dimz=', 'number=', 'prec=']
try:
    opts, args = getopt.getopt(sys.argv[1:], short_opts, long_opts)
except getopt.GetoptError:
    usage()
    sys.exit(2)

if len(args) != 0:
    filename = args[0]
    f = open(args[0], 'wt')
else:
    f = sys.stdout
    
# specify default and overwrite with new config
conf = { '--dimf': '2', '--dimz': '6,15', \
         '--number': '1000', '--prec': '3' }
         
conf.update(dict(opts))

if conf['--dimz'] != '0':
    cards = map(lambda x: int(x), conf['--dimz'].split(","))
else:
    cards = []

# print the number of discrete space dimensions
f.write("%d\n" % (len(cards)))

# print the discrete space dimensions' cardinalities
for c in cards:
    f.write("%d " % (c))
f.write("\n")

# write the number of continuous space dimensions
f.write("%d\n" % (int(conf['--dimf'])))

# print the number of points
f.write("%d\n" % (int(conf['--number'])))
format = "%%.%df " % (int(conf['--prec']))
# print n random points
for n in range(0, int(conf['--number'])):
    for c in cards:
        f.write("%d " % (int(random.random() * (c - 1) + 1)))
    for i in range(0, int(conf['--dimf'])):
        f.write(format % (random.random()))
    f.write("\n") 
