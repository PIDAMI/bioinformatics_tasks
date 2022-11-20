import numpy as np
import math
inf = math.inf
from profile_alingment import argmax


class PairwiseAligner(object):
    def __init__(self, **kwargs) -> None:
        self.delta = kwargs["delta"]
        self.eps = kwargs["eps"]
        self.tau = kwargs["tau"]
        self.q = kwargs["q"]
        self.p = kwargs["p"]

        self.s1 = kwargs["s1"]
        self.s2 = kwargs["s2"]
        self.n = len(self.s1) - 1
        self.m = len(self.s2) - 1
        self.res = [[[0,0,0] for j in range(self.n + 1)] for i in range(self.m + 1)]
        self.res[0][0][0] = 1
    
    def fit(self):
        for letter_num in range(1, self.n + 1):
            self.process_I(0, letter_num) # at 0th row no other states are valid
        for state_num in range(1, self.m + 1):
            self.process_D(state_num, 0)

        for i in range(1, self.m + 1):
            for j in range(1, self.n + 1):
                self.process_M(i, j)
                self.process_D(i, j)
                self.process_I(i, j)

    def process_M(self, state_num, letter_num):
        pass  
    def process_I(self, state_num, letter_num):
        pass  
    def process_D(self, state_num, letter_num): 
        pass 

    

## too lazy to rename states, I should be X and D should be named Y to match pairwise alingment notations
class ViterbiPairwiseAligner(PairwiseAligner):

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)
        self.pointers = [[[0,0,0] for j in range(self.n + 1)] for i in range(self.m + 1)]

        
    def fit(self):
        super().fit()
        score_index, score = argmax(self.res[-1][-1])
        self.score = score * self.tau
        self.score_index = score_index
        self.alingment = self.traceback_align()


    def traceback_align(self):
        
        aling1 = ''
        aling2 = ''
        letter_ind = self.n
        state_ind = self.m 
        state = self.score_index # M, I or D
        new_letters = self.get_letter(self.s1, self.s2, state, state_ind, letter_ind)
        aling1 += new_letters[0]
        aling2 += new_letters[1]
        
        while letter_ind > 0:
            state_ind, letter_ind,state = self.trace_one_step(state_ind, letter_ind, state)
            if letter_ind == 0:
                break
            new_letters = self.get_letter(self.s1, self.s2, state, state_ind, letter_ind)
            aling1 += new_letters[0]
            aling2 += new_letters[1]
        return aling1[::-1], aling2[::-1]

    
    def get_letter(self, seq1, seq2, state, state_num, letter_num):
        if state == 0: # match
            return seq1[letter_num], seq2[state_num] 
        elif state == 1: #insert
            return seq1[letter_num], "-"
        else: # del
            return "-",seq2[state_num]
            

    def trace_one_step(self, state_ind, letter_ind, state):
        res = [0, 0, self.pointers[state_ind][letter_ind][state]]
        if state == 0:
            res[0] = state_ind - 1
            res[1] = letter_ind - 1
        elif state == 1:
            res[0] = state_ind
            res[1] = letter_ind - 1
        else:
           res[0] = state_ind - 1
           res[1] = letter_ind
        return res


    def process_M(self, state_num, letter_num):  

        MM = (1 - 2 * self.delta - self.tau) * self.res[state_num-1][letter_num-1][0]
        IM = (1 - self.eps - self.tau) * self.res[state_num-1][letter_num-1][1]
        DM = (1 - self.eps - self.tau) * self.res[state_num-1][letter_num-1][2]
        ind,max = argmax([MM,IM,DM])    
        self.pointers[state_num][letter_num][0] = ind

        self.res[state_num][letter_num][0] = max * self.p[f'{self.s1[letter_num]}{self.s2[state_num]}']

    def process_I(self, state_num, letter_num):
       
        if state_num == 0 and letter_num == 1:
            ind = 0
            max = self.delta * self.res[state_num][letter_num-1][0]
        elif state_num  == 0:
            ind = 1
            max = self.eps * self.res[state_num][letter_num-1][1]
        else:
            MI = self.delta * self.res[state_num][letter_num-1][0]
            II = self.eps * self.res[state_num][letter_num-1][1]
            ind,max = argmax([MI,II])    
        self.pointers[state_num][letter_num][1] = ind
        self.res[state_num][letter_num][1] = max * self.q[f"{self.s1[letter_num]}"]

    def process_D(self, state_num, letter_num):
        if state_num == 1 and letter_num == 0:
            ind = 0
            max = self.delta * self.res[state_num-1][letter_num][0]
        elif letter_num  == 0:
            ind = 2
            max = self.eps * self.res[state_num-1][letter_num][2] 
        else:
            MD = self.delta * self.res[state_num-1][letter_num][0]
            DD = self.eps * self.res[state_num-1][letter_num][2]
            ind,max = argmax([MD,DD])    
        self.pointers[state_num][letter_num][2] = ind
        self.res[state_num][letter_num][2] = max * self.q[f"{self.s2[state_num]}"]

class ForwardPairwiseAligner(PairwiseAligner):

    def __init__(self, **kwargs) -> None:
        super().__init__(**kwargs)

        
    def fit(self):
        super().fit()
        score = sum(self.res[-1][-1])
        self.score = score * self.tau 


    def process_M(self, state_num, letter_num):  

        MM = (1 - 2 * self.delta - self.tau) * self.res[state_num-1][letter_num-1][0]
        IM = (1 - self.eps - self.tau) * self.res[state_num-1][letter_num-1][1]
        DM = (1 - self.eps - self.tau) * self.res[state_num-1][letter_num-1][2]
        
        self.res[state_num][letter_num][0] = sum([MM,IM,DM]) * self.p[f'{self.s1[letter_num]}{self.s2[state_num]}']

    def process_I(self, state_num, letter_num):
       
        if state_num == 0 and letter_num == 1:
            max = self.delta * self.res[state_num][letter_num-1][0]
        elif state_num  == 0:
            max = self.eps * self.res[state_num][letter_num-1][1]
        else:
            MI = self.delta * self.res[state_num][letter_num-1][0]
            II = self.eps * self.res[state_num][letter_num-1][1]
            max = MI + II
        self.res[state_num][letter_num][1] = max * self.q[f"{self.s1[letter_num]}"]

    def process_D(self, state_num, letter_num):
        if state_num == 1 and letter_num == 0:
            max = self.delta * self.res[state_num-1][letter_num][0]
        elif letter_num  == 0:
            max = self.eps * self.res[state_num-1][letter_num][2] 
        else:
            MD = self.delta * self.res[state_num-1][letter_num][0]
            DD = self.eps * self.res[state_num-1][letter_num][2]
            max = MD + DD   
        self.res[state_num][letter_num][2] = max * self.q[f"{self.s2[state_num]}"]



# delta = 0.2
# tau = 0.1
# eps = 0.1
# q = {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25}
# p = dict()
# for l1 in ["A", "C", "G", 'T']:
#     for l2 in ["A", "C", "G", 'T']:
#         val = 0
#         if l1 == l2:
#             val = 0.5
#         elif l1 in ["C", "T"] and l2 in ["C", "T"]:
#             val = 0.05
#         elif l1 in ["A", "G"] and l2 in ["A", "G"]:
#             val = 0.05
#         elif l1 in ["A", "T"] and l2 in ["A", "T"]:
#             val = 0.3
#         elif l1 in ["G", "C"] and l2 in ["G", "C"]:
#             val = 0.3
#         elif l1 in ["G", "T"] and l2 in ["G", "T"]: 
#             val = 0.15
#         elif l1 in ["A", "C"] and l2 in ["A", "C"]:
#             val = 0.15              
        
#         p[f"{l1}{l2}"] = val






