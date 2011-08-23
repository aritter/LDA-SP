#!/usr/bin/python

import sys
from TopicDistributions import *

class LdaEvaluator:
    nTopics = 300
    def __init__(self, dataDir=None):
        self.a1TopicWord = TopicWordDistributions(dataDir + "/arg1-topic-word")
        self.a2TopicWord = TopicWordDistributions(dataDir + "/arg2-topic-word")
        self.predTopic   = PredTopicDistributions(dataDir + "/doc-topic")
        self.a1Cnt = ArgCnt(dataDir + "/a1cnt")
        self.a2Cnt = ArgCnt(dataDir + "/a2cnt")

    def ProbNotDistractor1(self, pred, arg):
        pArgGivenPred = 0.0
        for i in range(1,self.nTopics+1):
            pArgGivenPred += self.a1TopicWord.ProbWordGivenTopic(arg, str(i)) * self.predTopic.ProbTopicGivenPred(str(i), pred)
        pArg = self.a1Cnt.GetProb(arg)
        pNotDistractor  = pArgGivenPred  * 0.5 / (0.5*pArg  + 0.5*pArgGivenPred)
        return pNotDistractor

    def ProbNotDistractor2(self, pred, arg):
        pArgGivenPred = 0.0
        for i in range(1,self.nTopics+1):
            pArgGivenPred += self.a2TopicWord.ProbWordGivenTopic(arg, str(i)) * self.predTopic.ProbTopicGivenPred(str(i), pred)
        pArg = self.a2Cnt.GetProb(arg)
        pNotDistractor  = pArgGivenPred  * 0.5 / (0.5*pArg  + 0.5*pArgGivenPred)
        return pNotDistractor

if __name__ == "__main__":
    ldaEval = LdaEvaluator('lda_model')
    for line in open(sys.argv[1]):
        line = line.rstrip('\n')
        (pred, arg, label) = line.split('\t')
        if int(sys.argv[2]) == 1:
            score = ldaEval.ProbNotDistractor1(pred, arg)
        elif int(sys.argv[2]) == 2:
            score = ldaEval.ProbNotDistractor2(pred, arg)
        print "%s\t%s\t%s" % (pred, arg, score)
