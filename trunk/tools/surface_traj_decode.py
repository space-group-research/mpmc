#!/usr/bin/python

import sys 

# KEY:
#
# + Newlines are signified via semicolons.
# + Replace "," with " " (CSV-like format)
# + "@" --> "ATOM"
# + "#" --> "REMARK step="
# + "Rn" --> "REMARK node="
# + "*" --> "0.000"
# + "!" --> "ENDMDL"
#

def decode(trajfile):

    # Open file, read into a string, then close
    traj = open(trajfile)
    print "Reading source file: " + trajfile
    tstr = traj.read()	# Read file (1 VERY long line) into a string for parsing
    traj.close()

    # Create output file
    of = "surface_trajectory.pqr"
    print "Writing output to file: " + of
    of = open(of, "w")

    # Generate newlines
    tlines = tstr.split(";")	# tlines is a list, not a string

    i = 0
    for line in tlines:
	space = line.replace(",", " ")
	zers  = space.replace("*", "0.000")
	end   = zers.replace("!", "ENDMDL")
	if end == "#":
	    i += 1
	    num = end.replace("#", str(i))
	    of.write("REMARK step=%s\n" % num)
	elif end[0] == "@":
	    atom  = end.replace("@", "ATOM")
	    burst = atom.split()
	    of.write("%s  %5d %-4.45s %-3.3s %-1.1s %4d   %8.3f%8.3f%8.3f\n" % (burst[0],int(burst[1]),burst[2],burst[3],burst[4],int(burst[5]),float(burst[6]),float(burst[7]),float(burst[8])))
	else:
	    of.write("%s\n" % end)

    # Clean up
    of.close()


if __name__ == "__main__":

    # Check args
    if (len(sys.argv) == 2):
	trajfile = sys.argv[1]
    else:
	print "usage: ./surface_traj_inflate.py [FILE]"
	sys.exit(1)

    # Decode file
    decode(trajfile)

