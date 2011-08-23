#!/usr/bin/python

import math
import sys
from TopicDistributions import *
from optparse import OptionParser

class LdaInferenceEvaluator:
    nTopics = 300
    def __init__(self, dataDir='lda_model'):
        self.a1TopicWord = TopicWordDistributions(dataDir + "/arg1-topic-word")
        self.a2TopicWord = TopicWordDistributions(dataDir + "/arg2-topic-word")
        self.predTopic   = PredTopicDistributions(dataDir + "/doc-topic")

    def ProbWord1(self, pred, arg):
        probSum = 0.0
        for i in range(1,self.nTopics+1):
            probSum += self.a1TopicWord.ProbWordGivenTopic(arg, str(i)) * self.predTopic.ProbTopicGivenPred(str(i), pred)
        return probSum

    def ProbWord2(self, pred, arg):
        probSum = 0.0
        for i in range(1,self.nTopics+1):
            sum += self.a2TopicWord.ProbWordGivenTopic(arg, str(i)) * self.predTopic.ProbTopicGivenPred(str(i), pred)
        return probSum

    def Evaluate(self, antecedent, consequent):
        #Only want to compute these once
        aProb1 = self.ProbWord1(antecedent[1], antecedent[0])
        aProb2 = self.ProbWord1(antecedent[1], antecedent[2])
        cProb1 = self.ProbWord1(consequent[1], consequent[0])
        cProb2 = self.ProbWord1(consequent[1], consequent[2])

        #TODO: do we need to use logs here?

        #Sum over all topic pairs
        probSum = 0.0
        for i in range(1,self.nTopics+1):
            for j in range(1,self.nTopics+1):
                #Probability that both antecedent and consequent have topic i for arg1 and topic j for arg2
                probSum += (self.a1TopicWord.ProbWordGivenTopic(antecedent[0], str(i)) * self.predTopic.ProbTopicGivenPred(str(i), antecedent[1]) *
                            self.a2TopicWord.ProbWordGivenTopic(antecedent[2], str(i)) * self.predTopic.ProbTopicGivenPred(str(i), antecedent[1]) *
                            self.a1TopicWord.ProbWordGivenTopic(consequent[0], str(i)) * self.predTopic.ProbTopicGivenPred(str(i), consequent[1]) *
                            self.a2TopicWord.ProbWordGivenTopic(consequent[2], str(i)) * self.predTopic.ProbTopicGivenPred(str(i), consequent[1]))
        return probSum

if __name__ == '__main__':
    ldaie = LdaInferenceEvaluator('lda_model')
    inferences = []
    for line in open(sys.argv[1]):
        line = line.rstrip('\n')
        line = line.rstrip('\t')
        (aArg1, aPred, aArg2, cArg1, cPred, cArg2, label) = line.split('\t')
        p = ldaie.Evaluate((aArg1, aPred, aArg2), (cArg1, cPred, cArg2))
        print "\t".join([aArg1, aPred, aArg2, cArg1, cPred, cArg2, str(math.log(p))])
