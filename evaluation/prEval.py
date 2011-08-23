#!/usr/bin/python

import sys

def PrEval(labelFiles, predFiles, threshold):
    tp = 0
    fp = 0
    tn = 0
    fn = 0
    for i in range(len(labelFiles)):
        predIn = open(predFiles[i])
        for line in open(labelFiles[i]):
            label = line.rstrip('\n').split('\t')[-1]
            pred  = float(predIn.readline().rstrip('\n').split('\t')[-1])

            if label == '+':
                if pred > threshold:
                    tp += 1
                else:
                    fn += 1
            elif label == '-':
                if pred > threshold:
                    fp += 1
                else:
                    tn += 1
            else:
                raise Exception('Improper labels')
        predIn.close()

    p = float(tp) / float(tp + fp)
    r = float(tp) / float(tp + fn)
    f = 2*p*r / (p + r)
    return (p,r,f)

if __name__ == "__main__":
    thresh = float(sys.argv[1])
    labelFiles = [x.split(':')[0] for x in sys.argv[2:]]
    predFiles  = [x.split(':')[1] for x in sys.argv[2:]]

    (p,r,f) = PrEval(labelFiles, predFiles, thresh)

    print "P:%s\nR:%s\nF:%s" % (p,r,f)
