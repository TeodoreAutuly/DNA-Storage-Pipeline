#### Libraries ####
import numpy as np
import os
import random

#### Functions ###
# Defines the class containing the model + the function to load the model
class Str_Model:
    def __init__(self, path):
	# Useful variables for the simulator
        self.subList = np.array([["A2C", "A2G", "A2T"],
                                 ["C2A", "C2G", "C2T"],
                                 ["G2A", "G2C", "G2T"],
                                 ["T2A", "T2C", "T2G"]])
        self.Yi = np.array(["M", "D", "I", "S"])

        # Load value of K
        with open(path + "/Kval.txt", "r") as file:
            self.Lkmer = int(file.readline().strip())

        # Load average edit probabilities in the beginning of the sequence
        self.BegProbEdit = np.array([float(x) for x in open(path + "/BegErrByPosAvg.txt").readline().strip().split()])
        self.BegProbIns = np.array([float(x) for x in open(path + "/BegInsByPosAvg.txt").readline().strip().split()])
        BegProbDel = np.array([float(x) for x in open(path + "/BegDelByPosAvg.txt").readline().strip().split()])
        BegProbSubst = np.array([float(x) for x in open(path + "/BegMisAvg.txt").readline().strip().split()])

        # Convert all these probability distributions into cumulative distribution functions
        if (BegProbDel.sum() + BegProbSubst.sum()) == 0:
            self.rangeDelMaxBeg = 0
            self.rangeSubstMaxBeg = 1
        else:
            sumBeg = BegProbDel + BegProbSubst
            self.rangeDelMaxBeg = BegProbDel / sumBeg
            self.rangeSubstMaxBeg = self.rangeDelMaxBeg + (BegProbSubst / sumBeg)

        # Load average edit probabilities in the end of the sequence
        self.EndProbEdit = np.array([float(x) for x in open(path + "/EndErrByPosAvg.txt").readline().strip().split()])
        self.EndProbIns = np.array([float(x) for x in open(path + "/EndInsByPosAvg.txt").readline().strip().split()])
        EndProbDel = np.array([float(x) for x in open(path + "/EndDelByPosAvg.txt").readline().strip().split()])
        EndProbSubst = np.array([float(x) for x in open(path + "/EndMisAvg.txt").readline().strip().split()])

        # Convert all these probability distributions into cumulative distribution functions
        if (EndProbDel.sum() + EndProbSubst.sum()) == 0:
            self.rangeDelMaxEnd = 0
            self.rangeSubstMaxEnd = 1
        else:
            sumEnd = EndProbDel + EndProbSubst
            self.rangeDelMaxEnd = EndProbDel / sumEnd
            self.rangeSubstMaxEnd = self.rangeDelMaxEnd + (EndProbSubst / sumEnd)

	# Load average edit probabilities per base A,C,G,T
        self.MidProbEditA = np.array([float(x) for x in open(path + "/AErrAvg.txt").readline().strip().split()])
        self.MidProbInsA = np.array([float(x) for x in open(path + "/AInsAvg.txt").readline().strip().split()])
        MidProbDelA = np.array([float(x) for x in open(path + "/ADelAvg.txt").readline().strip().split()])
        MidProbSubstA = np.array([float(x) for x in open(path + "/AMisAvg.txt").readline().strip().split()])

        self.MidProbEditC = np.array([float(x) for x in open(path + "/CErrAvg.txt").readline().strip().split()])
        self.MidProbInsC = np.array([float(x) for x in open(path + "/CInsAvg.txt").readline().strip().split()])
        MidProbDelC = np.array([float(x) for x in open(path + "/CDelAvg.txt").readline().strip().split()])
        MidProbSubstC = np.array([float(x) for x in open(path + "/CMisAvg.txt").readline().strip().split()])

        self.MidProbEditG = np.array([float(x) for x in open(path + "/GErrAvg.txt").readline().strip().split()])
        self.MidProbInsG = np.array([float(x) for x in open(path + "/GInsAvg.txt").readline().strip().split()])
        MidProbDelG = np.array([float(x) for x in open(path + "/GDelAvg.txt").readline().strip().split()])
        MidProbSubstG = np.array([float(x) for x in open(path + "/GMisAvg.txt").readline().strip().split()])

        self.MidProbEditT = np.array([float(x) for x in open(path + "/TErrAvg.txt").readline().strip().split()])
        self.MidProbInsT = np.array([float(x) for x in open(path + "/TInsAvg.txt").readline().strip().split()])
        MidProbDelT = np.array([float(x) for x in open(path + "/TDelAvg.txt").readline().strip().split()])
        MidProbSubstT = np.array([float(x) for x in open(path + "/TMisAvg.txt").readline().strip().split()])

        # Convert all these probability distributions into cumulative distribution functions, per base
        # A
        if (MidProbDelA.sum() + MidProbSubstA.sum()) == 0:
            self.rangeDelMaxMidA = 0
            self.rangeSubstMaxMidA = 1
        else:
            sumMidA = MidProbDelA + MidProbSubstA
            self.rangeDelMaxMidA = MidProbDelA / sumMidA
            self.rangeSubstMaxMidA = self.rangeDelMaxMidA + (MidProbSubstA / sumMidA)

        # C
        if (MidProbDelC.sum() + MidProbSubstC.sum()) == 0:
            self.rangeDelMaxMidC = 0
            self.rangeSubstMaxMidC = 1
        else:
            sumMidC = MidProbDelC + MidProbSubstC
            self.rangeDelMaxMidC = MidProbDelC / sumMidC
            self.rangeSubstMaxMidC = self.rangeDelMaxMidC + (MidProbSubstC / sumMidC)

        # G
        if (MidProbDelG.sum() + MidProbSubstG.sum()) == 0:
            self.rangeDelMaxMidG = 0
            self.rangeSubstMaxMidG = 1
        else:
            sumMidG = MidProbDelG + MidProbSubstG
            self.rangeDelMaxMidG = MidProbDelG / sumMidG
            self.rangeSubstMaxMidG = self.rangeDelMaxMidG + (MidProbSubstG / sumMidG)

        # T
        if (MidProbDelT.sum() + MidProbSubstT.sum()) == 0:
            self.rangeDelMaxMidT = 0
            self.rangeSubstMaxMidT = 1
        else:
            sumMidT = MidProbDelT + MidProbSubstT
            self.rangeDelMaxMidT = MidProbDelT / sumMidT
            self.rangeSubstMaxMidT = self.rangeDelMaxMidT + (MidProbSubstT / sumMidT)

	# Load average insertion length, in the beginning, per base A,C,G,T, and in the end
        probInsLenBeg = np.array([float(x) for x in open(path + "/insLenBegRates.txt").readline().strip().split()[:-1]])
        probInsLenMidA = np.array([float(x) for x in open(path + "/AinsLenMidRates.txt").readline().strip().split()[:-1]])
        probInsLenMidC = np.array([float(x) for x in open(path + "/CinsLenMidRates.txt").readline().strip().split()[:-1]])
        probInsLenMidG = np.array([float(x) for x in open(path + "/GinsLenMidRates.txt").readline().strip().split()[:-1]])
        probInsLenMidT = np.array([float(x) for x in open(path + "/TinsLenMidRates.txt").readline().strip().split()[:-1]])
        probInsLenEnd = np.array([float(x) for x in open(path + "/insLenEndRates.txt").readline().strip().split()[:-1]])

        # Converts all these probabilities in cumulative distribution functions
        #self.rangeLenInsBeg = np.cumsum(probInsLenBeg) / np.sum(probInsLenBeg)
        #self.rangeLenInsMidA = np.cumsum(probInsLenMidA) / np.sum(probInsLenMidA)
        #self.rangeLenInsMidC = np.cumsum(probInsLenMidC) / np.sum(probInsLenMidC)
        #self.rangeLenInsMidG = np.cumsum(probInsLenMidG) / np.sum(probInsLenMidG)
        #self.rangeLenInsMidT = np.cumsum(probInsLenMidT) / np.sum(probInsLenMidT)
        #self.rangeLenInsEnd = np.cumsum(probInsLenEnd) / np.sum(probInsLenEnd)

        sum_probInsLenBeg = np.sum(probInsLenBeg)
        if sum_probInsLenBeg == 0:
            self.rangeLenInsBeg = 1
        else:
            self.rangeLenInsBeg = np.cumsum(probInsLenBeg) / sum_probInsLenBeg

        sum_probInsLenMidA = np.sum(probInsLenMidA)
        if sum_probInsLenMidA == 0:
            self.rangeLenInsMidA = [1,0]
        else:
            self.rangeLenInsMidA = np.cumsum(probInsLenMidA) / sum_probInsLenMidA

        sum_probInsLenMidC = np.sum(probInsLenMidC)
        if sum_probInsLenMidC == 0:
            self.rangeLenInsMidC = [1,0]
        else:
            self.rangeLenInsMidC = np.cumsum(probInsLenMidC) / sum_probInsLenMidC

        sum_probInsLenMidG = np.sum(probInsLenMidG)
        if sum_probInsLenMidG == 0:
            self.rangeLenInsMidG = [1,0]
        else:
            self.rangeLenInsMidG = np.cumsum(probInsLenMidG) / sum_probInsLenMidG

        sum_probInsLenMidT = np.sum(probInsLenMidT)
        if sum_probInsLenMidT == 0:
            self.rangeLenInsMidT = [1,0]
        else:
            self.rangeLenInsMidT = np.cumsum(probInsLenMidT) / sum_probInsLenMidT

        sum_probInsLenEnd = np.sum(probInsLenEnd)
        if sum_probInsLenEnd == 0:
            self.rangeLenInsEnd = [1,0]
        else:
            self.rangeLenInsEnd = np.cumsum(probInsLenEnd) / sum_probInsLenEnd

	# Construct transition probability matrix
        self.transProb = np.zeros((4, 3))
        self.rangeTransProb = np.zeros((4, 3))
        for i in range(len(self.subList)):
            # transProb
            for j in range(len(self.subList[i])):
                self.transProb[i, j] = float(open(path + f"/{self.subList[i][j]}_Avg.txt").readline().strip())
            self.transProb[i, :] = self.transProb[i, :] / np.sum(self.transProb[i, :])

            # rangeTransProb
            self.rangeTransProb[i, 0] = self.transProb[i, 0]
            for j in range(1, len(self.subList[i])):
                self.rangeTransProb[i, j] = self.rangeTransProb[i, j - 1] + self.transProb[i, j]

	#Probability to observe event Ei on kmer i, knowing that event Ei-1 occured
        #index: 1-> prev Match, 2-> prev Del, 3-> prev Ins, 4-> prev Subst,
        #Yi-1= [Match[Yi=I,D,S,E,M],Del[Yi=I,D,S,E,M],Ins[Yi=I,D,S,E,M],Sub[Yi=I,D,S,E,M]]
        self.EditYiKmerprevYi = []
        self.mapKmerPrevYi = [{} for _ in range(4)]  # Array of dictionaries

        for e1 in range(4):
            if os.path.isfile(path + f"/KmerYi_prevYi{self.Yi[e1]}_RatesAvg.txt"):
                # kmer Ins Del Subst Err Match
                with open(path + f"/KmerYi_prevYi{self.Yi[e1]}_RatesAvg.txt", "r") as file:
                     edit_yi_kmer_prev_yi = np.genfromtxt(file, dtype=None, names=None, encoding=None)
                nKmer = edit_yi_kmer_prev_yi.shape[0]  # Number of entries

                # Create a dictionary (kmer,index) to reach the index (in KmerProb) containing probabilities (edit Rates table)
                for i in range(nKmer):
                    self.mapKmerPrevYi[e1][edit_yi_kmer_prev_yi[i][0]] = i
                self.EditYiKmerprevYi.append(edit_yi_kmer_prev_yi)
            else:
                self.EditYiKmerprevYi.append([])
                self.mapKmerPrevYi[e1] = {}


        # Probability of insertion length, depending on the K-mer and on previous error event
        self.kmerInsLenProbPrevYi = []
        self.mapKmerInsLenPrevYi = [{} for _ in range(4)]

        for e1 in range(4):
               if os.path.isfile(path + f"/KmerInsLen_prevYi{self.Yi[e1]}_RatesAvg2.txt"):
                     with open(path + f"/KmerInsLen_prevYi{self.Yi[e1]}_RatesAvg2.txt", "r") as file:
                          kmer_ins_len_prob_tmp = np.genfromtxt(file, dtype=None, names=None, encoding=None)
                     if(kmer_ins_len_prob_tmp.ndim == 1):
                          n_kmer = 1
                     else:
                          n_kmer = kmer_ins_len_prob_tmp.shape[0]
                     kmer_ins_range_prob_tmp = []

                     for i in range(n_kmer):
                          if(n_kmer==1):
                              floats = list(map(float, kmer_ins_len_prob_tmp[1].split(",")))
                              self.mapKmerInsLenPrevYi[e1][kmer_ins_len_prob_tmp[0]] = i
                          else:
                              floats = list(map(float, kmer_ins_len_prob_tmp[i][1].split(",")))
                              self.mapKmerInsLenPrevYi[e1][kmer_ins_len_prob_tmp[i][0]] = i
                          tmp_line_cumsum = np.cumsum(floats)
                          kmer_ins_range_prob_tmp.append(tmp_line_cumsum)
                     self.kmerInsLenProbPrevYi.append(kmer_ins_range_prob_tmp)

               else:
                     self.mapKmerInsLenPrevYi[e1] = {}
                     self.kmerInsLenProbPrevYi[e1] = []

    def channel_simulator(self, Xdna, nbrOut):
        # Initialization
        nRef = len(Xdna)
        SimSeq = []

        # In the simulator, the following convention is used for events
        # 1 : Match
        # 2 : Del
        # 3 : Ins
        # 4 : Sub

        for iSeq in range(0, nbrOut):
            ####### Builds sequence of events ######
            Events = np.zeros(nRef, dtype=int)       # Events
            Ins_length = np.zeros(nRef, dtype=int)   # When insertion, insertion length
            Sub_event = [] # When substitution (or insertion/substitution or insertion/match), replace by which base

            ### Generates successive error events
            ### First position
            E, L, Bval = self.sample_event(self.BegProbEdit, self.BegProbIns, self.rangeDelMaxBeg,
                                            self.rangeSubstMaxBeg, self.rangeLenInsBeg, self.rangeTransProb,
                                            self.subList, ''.join(Xdna[0]))
            Events[0] = E
            Ins_length[0] = L
            Sub_event.append(Bval)

            ### Less than k bases in the nanopore => Probability only depends on the current base
            for iRef in range(1, self.Lkmer):
                current_base = ''.join(Xdna[iRef])
                if current_base == "A":
                    Events[iRef], Ins_length[iRef], Bval = self.sample_event(self.MidProbEditA, self.MidProbInsA,
                                                                              self.rangeDelMaxMidA, self.rangeSubstMaxMidA,
                                                                              self.rangeLenInsMidA, self.rangeTransProb,
                                                                              self.subList, ''.join(Xdna[iRef]))
                    Sub_event.append(Bval[0])
                elif current_base == "C":
                    Events[iRef], Ins_length[iRef], Bval = self.sample_event(self.MidProbEditC, self.MidProbInsC,
                                                                              self.rangeDelMaxMidC, self.rangeSubstMaxMidC,
                                                                              self.rangeLenInsMidC, self.rangeTransProb,
                                                                              self.subList, ''.join(Xdna[iRef]))
                    Sub_event.append(Bval[0])
                elif current_base == "T":
                    Events[iRef], Ins_length[iRef], Bval = self.sample_event(self.MidProbEditT, self.MidProbInsT,
                                                                              self.rangeDelMaxMidT, self.rangeSubstMaxMidT,
                                                                              self.rangeLenInsMidT, self.rangeTransProb,
                                                                              self.subList, ''.join(Xdna[iRef]))
                    Sub_event.append(Bval[0])
                elif current_base == "G":
                    Events[iRef], Ins_length[iRef], Bval = self.sample_event(self.MidProbEditG, self.MidProbInsG,
                                                                              self.rangeDelMaxMidG, self.rangeSubstMaxMidG,
                                                                              self.rangeLenInsMidG, self.rangeTransProb,
                                                                              self.subList, ''.join(Xdna[iRef]))
                    Sub_event.append(Bval[0])
                else:
                    print(Xdna[iRef])
                    raise ValueError("Unknown letter in the input")


            ### more than k bases passed through the nanopore
            for iRef in range(self.Lkmer, nRef-1):
                # Current K-mer
                kmer = ''.join(Xdna[iRef-self.Lkmer+1:iRef])

                # Reads probabilities for the current K-mer
                if Events[iRef-1] in self.mapKmerPrevYi and kmer in self.mapKmerPrevYi[Events[iRef-1]]:
                    indexP = self.mapKmerPrevYi[Events[iRef-1]][kmer]
                else:  # unknown kmer -> use average values
                    indexP = self.mapKmerPrevYi[Events[iRef-1]-1]["AVG"]

                probInsK = self.EditYiKmerprevYi[Events[iRef-1]-1][indexP][1]
                probDelK = self.EditYiKmerprevYi[Events[iRef-1]-1][indexP][2]
                probSubstK = self.EditYiKmerprevYi[Events[iRef-1]-1][indexP][3]
                probMatchK = self.EditYiKmerprevYi[Events[iRef-1]-1][indexP][5]
                probEditK = probDelK + probSubstK

                # Read insertion length probability for the current K-mer
                if Events[iRef-1] in self.mapKmerInsLenPrevYi and kmer in self.mapKmerInsLenPrevYi[Events[iRef-1]]:
                    indexPLen = self.mapKmerInsLenPrevYi[Events[iRef-1]-1][kmer]
                else:  # unknown -> use avg values (of all known kmers)
                    indexPLen = self.mapKmerInsLenPrevYi[Events[iRef-1]-1]["AVG"]

                rangeLenInsK = self.kmerInsLenProbPrevYi[Events[iRef-1]-1][indexPLen]

                # Sample error event
                rangeDelMaxK = probDelK / probEditK
                rangeSubstMaxK = rangeDelMaxK + (probSubstK / probEditK)
                Events[iRef], Ins_length[iRef], Bval = self.sample_event(probEditK, probInsK, rangeDelMaxK, rangeSubstMaxK,
                                                                     rangeLenInsK, self.rangeTransProb, self.subList,
                                                                     ''.join(Xdna[iRef]))
                Sub_event.append(Bval)

            ### Last position
            E, L, Bval = self.sample_event(self.EndProbEdit, self.EndProbIns, self.rangeDelMaxEnd,
                                        self.rangeSubstMaxEnd, self.rangeLenInsEnd, self.rangeTransProb,
                                        self.subList, ''.join(Xdna[nRef-1]))
            Events[nRef-1] = E
            Ins_length[nRef-1] = L
            Sub_event.append(Bval)

            ###### Generates output sequence ######
            CurrentSeq = []
            for iRef in range(nRef):
                if Events[iRef] == 1:  # Match
                    CurrentSeq.append(Xdna[iRef])
                elif Events[iRef] == 3:  # Insertion
                    CurrentSeq.append(Sub_event[iRef])
                    for _ in range(Ins_length[iRef]):
                        CurrentSeq.append(self.randomBase())
                elif Events[iRef] == 4:  # Substitution
                    CurrentSeq.append(Sub_event[iRef])
                elif Events[iRef] < 1 or Events[iRef] > 4:
                    raise ValueError("Got an unexpected event")
            SimSeq.append(CurrentSeq)

        return SimSeq


    def sample_event(self,probEdit, probIns, rangeDel, rangeSub, rangeLenIns, rangeTransProb, SubList, currentB):
        # Initialization
        r1 = random.random()
        r2 = random.random()
        r3 = random.random()

        # Error event
        if r1 <= probEdit:
            if r2 <= rangeDel:  # Deletion
                E = 2
                L = 0
                B = currentB
            elif r2 <= rangeSub:  # Substitution
                # Choose replacing base
                B = self.substB(rangeTransProb, currentB, SubList)
                if r3 <= probIns:  # Substitution + Insertion
                    E = 3
                    L = self.insertB(rangeLenIns)
                else:
                    E = 4
                    L = 0
            else:
                raise ValueError("Input range does not sum to 1\n")
        else:  # Match
            B = currentB
            if r3 <= probIns:  # Match + Insertion
                E = 3
                L = self.insertB(rangeLenIns)
            else:
                E = 1
                L = 0

        return E, L, B



    def substB(self, rangeTransProb, currBase, subList):
        # Initialization
        r = random.random()
        base = ""

        # Sample a base
        if currBase == "A":
            for j in range(len(rangeTransProb[0])):
                if r <= rangeTransProb[0][j]:
                    base = subList[0][j].split("2")[1]
                    break
        elif currBase == "C":
            for j in range(len(rangeTransProb[1])):
                if r <= rangeTransProb[1][j]:
                    base = subList[1][j].split("2")[1]
                    break
        elif currBase == "G":
            for j in range(len(rangeTransProb[2])):
                if r <= rangeTransProb[2][j]:
                    base = subList[2][j].split("2")[1]
                    break
        elif currBase == "T":
            for j in range(len(rangeTransProb[3])):
                if r <= rangeTransProb[3][j]:
                    base = subList[3][j].split("2")[1]
                    break
        else:
            print(currBase)
            raise ValueError("Unknown base in the function for base substitution\n")

        return base

    def insertB(self, rangeLenIns):
        # Initialization
        r3 = random.random()
        length = 1
        at = 0

        # Sample length of the insertion
        for j in range(1,len(rangeLenIns)):
            if r3 <= rangeLenIns[j]:
                length = j + 1
                at = 1
                break
        if at == 0:
            raise ValueError("Insertion probability does not sum to 1")

        return length

    def randomBase(self):
        base = random.choice(['A', 'C', 'G', 'T'])
        return base
