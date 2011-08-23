import re

def incHash2(hashtable, k1, k2, count):
    if not hashtable.has_key(k1):
        hashtable[k1] = {}
    if not hashtable[k1].has_key(k2):
        hashtable[k1][k2] = 0.0
    hashtable[k1][k2] += count

def incHash(hashtable, k1, count):
    if not hashtable.has_key(k1):
        hashtable[k1] = 0.0
    hashtable[k1] += count

class ArgCnt:
    alpha = 0.1 * 300.0
    def __init__(self, countFile):
        self.pArg = {}
        self.total = 0
        for line in open(countFile):
            line = line.rstrip('\n')
            (arg, count) = line.split('\t')
            count = int(count)
            self.pArg[arg] = count
            self.total += count

    def GetProb(self, arg):
        return float(self.pArg.get(arg,0) + self.alpha) / float(self.total + self.alpha * float(len(self.pArg.keys())))

class TopicWordDistributions:
    nTopics = 300
    eta = 0.1
    def __init__(self, topicFile):
        self.topicWordDistributions = {}
        self.topicTotals = {}
        self.vocab = {}

        for line in open(topicFile):
            line = line.rstrip('\n')
            topicWordData = line.split('\t')
            topic = topicWordData[0]
            self.topicWordDistributions[topic] = {}

            for i in range(1,len(topicWordData)):
                m = re.match(r'(.+):(\d+)', topicWordData[i])
                if not m:
                    #print topicWordData[i]
                    continue
                word = m.group(1)
                count = float(m.group(2))

                if not self.topicWordDistributions[topic].has_key(word):
                    self.topicWordDistributions[topic][word] = 0.0
                self.topicWordDistributions[topic][word] += count
                if not self.topicTotals.has_key(topic):
                    self.topicTotals[topic] = 0.0
                self.topicTotals[topic] += count
                self.vocab[word] = 1
        self.vocabSize = float(len(self.vocab.keys()))

    def ProbWordGivenTopic(self, word, topic):
        count = 0.0
        if self.topicWordDistributions[topic].has_key(word):
            count = self.topicWordDistributions[topic][word]
        return float(count + self.eta) / float(self.topicTotals[topic] + self.vocabSize * self.eta)

##########################################################################
# Note: this was copied from TextualInference.py (should have
# been in it's own file to begin with (maybe some refactoring is needed)
##########################################################################
class PredTopicDistributions:
    alpha = 0.1
    maxTopics = {}
    def __init__(self, inFile, nTopics=300):
        self.nTopics = nTopics 
        self.predTopicDistributions = {}
        self.predTotals = {}

        #Also maintain seperate distributions for arg1/arg2 (maybe this will improve performance?):
        self.arg1TopicDistributions = {}
        self.arg2TopicDistributions = {}
        self.arg1TopicTotals = {}
        self.arg2TopicTotals = {}

        for line in open(inFile):
            line = line.rstrip('\n')
            predTopicData = line.split('\t')
            pred = predTopicData[0]
            self.predTopicDistributions[pred] = {}
            self.arg1TopicDistributions[pred] = {}
            self.arg2TopicDistributions[pred] = {}

            for i in range(1,len(predTopicData)):
                m = re.match(r'(\d+)-(\d+):(\d+)', predTopicData[i])
                if not m:
                    continue
                t1 = m.group(1)
                t2 = m.group(2)
                count = float(m.group(3))

                #print "%s\t%s\t%f" % (t1, t2, count)

                #if self.highPMItopics == None or self.highPMItopics.has_key(t1):
                incHash2(self.predTopicDistributions, pred, t1, count)
                incHash2(self.predTopicDistributions, pred, t2, count)
                incHash(self.predTotals, pred, count)

                #if self.highPMItopics == None or self.highPMItopics.has_key(t2):
                incHash2(self.arg1TopicDistributions, pred, t1, count)
                incHash2(self.arg2TopicDistributions, pred, t2, count)
                incHash(self.predTotals, pred, count)

        #Find the topic with maximum probability
        self.maxCount = {}
        for pred in self.predTopicDistributions.keys():
            for topic in self.predTopicDistributions[pred].keys():
                if self.predTopicDistributions[pred][topic] > self.maxCount.get(pred,0):
                    self.maxTopics[pred] = topic
                    self.maxCount[pred] = self.predTopicDistributions[pred][topic]

    #Is the most probable topic for this predicate high PMI?
    def IsTopHighPMI(self, pred, arg):
        if arg == 1:
            argTopicPMI = self.arg1topicPMI
        if arg == 2:
            argTopicPMI = self.arg2topicPMI
        return argTopicPMI[self.maxTopics[pred]] > self.PMIthreshold

    def ProbTopicGivenPred(self, topic, pred):
        count = 0.0
        if self.predTopicDistributions[pred].has_key(topic):
            count = self.predTopicDistributions[pred][topic]

        #if self.highPMItopics != None and not self.highPMItopics.has_key(topic):
        #    return 0.0

        #print "%s\t%s\t%s" % (topic, pred, count)
        return float(count + self.alpha) / float(self.predTotals[pred] + float(self.nTopics) * self.alpha)

    def ProbTopicGivenPredArg(self, topic, pred, arg):
        count = 0.0
        if arg == 1:
            dist = self.arg1TopicDistributions
        elif arg == 2:
            dist = self.arg2TopicDistributions
            
        if dist[pred].has_key(topic):
            count = dist[pred][topic]
        return float(count + self.alpha) / float(self.predTotals[pred] / 2.0 + float(self.nTopics) * self.alpha)

    def WriteMaxTopics(self, outFile):
        fOut = open(outFile, 'w')
        for pred in self.maxTopics.keys():
            fOut.write("%s\t%s\n" % (pred, self.maxTopics[pred]))
