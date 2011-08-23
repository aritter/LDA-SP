#!/usr/bin/python

import sys
import prEval

## {{{ http://code.activestate.com/recipes/66472/ (r1)
def frange(start, end=None, inc=None):
    "A range function, that does accept float increments..."

    if end == None:
        end = start + 0.0
        start = 0.0

    if inc == None:
        inc = 1.0

    L = []
    while 1:
        next = start + len(L) * inc
        if inc > 0 and next >= end:
            break
        elif inc < 0 and next <= end:
            break
        L.append(next)
        
    return L
## end of http://code.activestate.com/recipes/66472/ }}}

if __name__ == "__main__":
    labelFiles = [x.split(':')[0] for x in sys.argv[1:]]
    predFiles  = [x.split(':')[1] for x in sys.argv[1:]]

    minScore = None
    maxScore = None
    for pf in predFiles:
        for line in open(pf):
            score = line.rstrip('\n').split('\t')[-1]
            if minScore == None or float(score) < minScore:
                minScore = float(score)
            if maxScore == None or float(score) > maxScore:
                maxScore = float(score)

    for threshold in frange(minScore, maxScore, (maxScore - minScore) / 10000):
        (p,r,f) = prEval.PrEval(labelFiles, predFiles, threshold)
        print "%s\t%s\t%s\t%s" % (threshold,p,r,f)
