#!/usr/bin/python

import sys
from optparse import OptionParser

parser = OptionParser()
parser.add_option("--dataDir", dest="dataDir", default=None,
                  help="directory containing tuples")
parser.add_option("--sampleFile", dest="sampleFile", default=None,
                  help="file containing samples")
parser.add_option("--nSamples", dest="nSamples", default=1,
                  help="number of samples")
parser.add_option("--tuplesOut", dest="tuplesOut", default=None,
                  help="output file for tagged tuples")
(options, args) = parser.parse_args()

print options.dataDir

##########################################################
# Read in vocabularies
##########################################################
a1vocab = {}
for line in open(options.dataDir + "/a1vocab"):
    line = line.rstrip('\n')
    (word, wid) = line.split('\t')
    a1vocab[wid] = word
a2vocab = {}
for line in open(options.dataDir + "/a2vocab"):
    line = line.rstrip('\n')
    (word, wid) = line.split('\t')
    a2vocab[wid] = word


##########################################################
# First pass: figure out how many samples total
##########################################################
totalSamples = 0
nPreds = 0

for line in open(options.sampleFile):
    if line[0:4] == "z2 =":
        totalSamples += 1
for line in open(options.dataDir + "/preds"):
    nPreds += 1

print totalSamples
print nPreds

##########################################################
# Second pass: compute distributions
##########################################################
docTopicCnt = {}
topicArg1Cnt = {}
topicArg2Cnt = {}

if options.tuplesOut:
    fTupOut = open(options.tuplesOut, 'w')

nProcessed = 0
z1 = None
z2 = None
for line in open(options.sampleFile):
    if line[0:4] == "z2 =":
        nProcessed += 1
    if totalSamples - nProcessed - 1 >= int(options.nSamples):
        # not ready to start taking samples yet
        continue

    #Which sample is this?
    line = line.rstrip('\n')
    if line[0:4] == "z1 =":
        z1 = [x.rstrip(' ').split(' ') for x in line[5:].split(' ; ')]
    elif line[0:4] == "z2 =":
        z2 = [x.rstrip(' ').split(' ') for x in line[5:].split(' ; ')]
        nProcessed += 1
    else:
        continue

    if z1 == None or z2 == None:
        continue

    ##########################################################
    # OK, we have samples for one doc read in now...
    ##########################################################

    #Read through arg1/arg2
    predsFile = open(options.dataDir + "/preds")
    a1File = open(options.dataDir + "/arg1")
    a2File = open(options.dataDir + "/arg2")
    for i in range(nPreds):
        pred = predsFile.readline().rstrip('\n')
        if not docTopicCnt.has_key(pred):
            docTopicCnt[pred] = {}

        arg1 = a1File.readline().rstrip('\n').split(' ')
        arg2 = a2File.readline().rstrip('\n').split(' ')
        for j in range(len(arg1)):
            if options.tuplesOut:
                fTupOut.write("%s\t%s\t%s\t%s\t%s\n" % (a1vocab[arg1[j]], z1[i][j], pred, a2vocab[arg2[j]], z2[i][j]))
            docTopicCnt[pred]["%s-%s" % (z1[i][j], z2[i][j])] = docTopicCnt[pred].get("%s-%s" % (z1[i][j], z2[i][j]), 0) + 1
            if not topicArg1Cnt.has_key(z1[i][j]):
                topicArg1Cnt[z1[i][j]] = {}
            topicArg1Cnt[z1[i][j]][arg1[j]] = topicArg1Cnt[z1[i][j]].get(arg1[j],0) + 1
            if not topicArg2Cnt.has_key(z2[i][j]):
                topicArg2Cnt[z2[i][j]] = {}
            topicArg2Cnt[z2[i][j]][arg2[j]] = topicArg2Cnt[z2[i][j]].get(arg2[j],0) + 1
        i += 1

    z1 = None
    z2 = None

##########################################################
# Finally, dump out distributions
##########################################################
def WriteDictSorted(inDict, outFile, mapper=lambda x: x, comparer=cmp):
    fOut = open(outFile, 'w')
    for key in sorted(inDict.keys(), lambda a,b: comparer(a,b)):
        fOut.write(str(key) + "\t" + "\t".join(["%s:%d" % (mapper(x),inDict[key][x]) for x in sorted(inDict[key].keys(), key=lambda k: inDict[key][k], reverse=True)]) + "\n")
    fOut.close()

WriteDictSorted(docTopicCnt,  "doc-topic")
WriteDictSorted(topicArg1Cnt, "arg1-topic-word", mapper=lambda x: a1vocab[x], comparer=lambda a,b: cmp(int(a), int(b)))
WriteDictSorted(topicArg2Cnt, "arg2-topic-word", mapper=lambda x: a2vocab[x], comparer=lambda a,b: cmp(int(a), int(b)))
