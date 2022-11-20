import numpy as np
import math
inf = math.inf

# return  max value as well as its position 
def argmax(vals):
    indexed_vals = list(enumerate(vals))
    return sorted(indexed_vals, key=lambda x: x[1])[-1]

ind_to_state = {0: "M", 1 : "I", 2 : "D"}



class Alinger(object):

    def __init__(self, **kwargs):
        self.trans_matrix = kwargs["trans_matrix"]
        self.gen_insert = kwargs["gen_insert"]
        self.gen_match = kwargs["gen_match"]
        self.seq = kwargs["seq"]
        self.q = kwargs["q"]
        self.n = len(self.seq) - 1 # seq size with respect to slack 1st letter
        self.m = len(self.trans_matrix) - 1  # number of Ms again with respect to slack var
        self.res = ([[[0,0,0] for i in range(self.n + 1)] for j in range(self.m + 1)])
        self.res[0][0] = [0,-inf,-inf]
        for i in range(1,self.n + 1):
            self.res[0][i] = [-inf,0,-inf]
        for j in range(1,self.m + 1):
            self.res[j][0] = [-inf,-inf,0]

    def process_M(self,state_num, letter_num):
        pass
    def process_D(self,state_num, letter_num):
        pass
    def process_I(self,state_num, letter_num):
        pass

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
        
    


class ViterbiProfileAligner(Alinger):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.pointers = [[[0,0,0] for i in range(self.n + 1)] for j in range(self.m + 1)]
        
    def fit(self):
        super().fit()
        score_index, score = argmax(self.res[-1][-1])
        self.score_index = score_index # 0,1 or 2 for M, I or D
        self.score = np.exp(score) * self.trans_matrix[-1][f"{ind_to_state[score_index]}M"] 
        self.alingment, self.states = self.traceback_align()
        
    
    def traceback_align(self):
        
        states = []
        aling = ''
        letter_ind = self.n
        state_ind = self.m 
        state = self.score_index # M, I or D
        states.insert(0,f"{ind_to_state[state]}{state_ind}")
        aling += self.get_letter(self.seq, state, letter_ind)
        while letter_ind > 0:
            state_ind, letter_ind,state = self.trace_one_step(state_ind, letter_ind, state)
            if letter_ind == 0:
                break
            states.insert(0,f"{ind_to_state[state]}{state_ind}")
            aling += self.get_letter(self.seq, state, letter_ind)
        return aling[::-1], states

    
    def get_letter(self, seq, state, letter_num):
        if state == 0: # match
            return seq[letter_num]
        elif state == 1: #insert
            return seq[letter_num]
        else: # delete
            return '-'

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
        if state_num == 1 and letter_num == 1:
            max = np.log(self.trans_matrix[state_num-1]["MM"]) + self.res[state_num-1][letter_num-1][0]
            self.pointers[state_num][letter_num][0] = 0
        elif state_num == 1:
            max = np.log(self.trans_matrix[state_num-1]["IM"]) + self.res[state_num-1][letter_num-1][1]
            self.pointers[state_num][letter_num][0] = 1

        elif letter_num == 1:
            max = np.log(self.trans_matrix[state_num-1]["DM"]) + self.res[state_num-1][letter_num-1][2]
            self.pointers[state_num][letter_num][0] = 2

        else:
            MM = np.log(self.trans_matrix[state_num-1]["MM"]) + self.res[state_num-1][letter_num-1][0]
            IM = np.log(self.trans_matrix[state_num-1]["IM"]) + self.res[state_num-1][letter_num-1][1]
            DM = np.log(self.trans_matrix[state_num-1]["DM"]) + self.res[state_num-1][letter_num-1][2]
            ind,max = argmax([MM,IM,DM])    
            self.pointers[state_num][letter_num][0] = ind

        letter = self.seq[letter_num]
        self.res[state_num][letter_num][0] = np.log(self.gen_match[state_num][letter] / self.q[letter]) + max

    def process_I(self, state_num, letter_num):
        if state_num == 0 and letter_num == 1:
            max = np.log(self.trans_matrix[state_num]["MI"]) + self.res[state_num][letter_num-1][0]
            self.pointers[state_num][letter_num][1] = 0
        elif letter_num == 1 and state_num != 0:
            max = np.log(self.trans_matrix[state_num]["DI"]) + self.res[state_num][letter_num-1][2]
            self.pointers[state_num][letter_num][1] = 2
        elif state_num == 0 and letter_num != 1:
            max = np.log(self.trans_matrix[state_num]["II"]) + self.res[state_num][letter_num-1][1]
            self.pointers[state_num][letter_num][1] = 1

        else:
            MI = np.log(self.trans_matrix[state_num]["MI"]) + self.res[state_num][letter_num-1][0]
            II = np.log(self.trans_matrix[state_num]["II"]) + self.res[state_num][letter_num-1][1]
            DI = np.log(self.trans_matrix[state_num]["DI"]) + self.res[state_num][letter_num-1][2]
            ind,max = argmax([MI,II,DI])    
            self.pointers[state_num][letter_num][1] = ind
        letter = self.seq[letter_num]
        self.res[state_num][letter_num][1] = np.log(self.gen_insert[state_num][letter] / self.q[letter]) + max

    def process_D(self, state_num, letter_num):
        if state_num == 1 and letter_num == 0:
            max =  np.log(self.trans_matrix[state_num-1]["MD"]) + self.res[state_num-1][letter_num][0]
            self.pointers[state_num][letter_num][2] = 0
        elif state_num == 1 and letter_num != 0:
            max = np.log(self.trans_matrix[state_num-1]["ID"]) + self.res[state_num-1][letter_num][1]
            self.pointers[state_num][letter_num][2] = 1
        elif state_num != 1 and letter_num == 0:
            max = np.log(self.trans_matrix[state_num-1]["DD"]) + self.res[state_num-1][letter_num][2]
            self.pointers[state_num][letter_num][2] = 2
        else:
            MD = np.log(self.trans_matrix[state_num-1]["MD"]) + self.res[state_num-1][letter_num][0]
            ID = np.log(self.trans_matrix[state_num-1]["ID"]) + self.res[state_num-1][letter_num][1]
            DD = np.log(self.trans_matrix[state_num-1]["DD"]) + self.res[state_num-1][letter_num][2]
            ind,max = argmax([MD,ID,DD])    
            self.pointers[state_num][letter_num][2] = ind


        self.res[state_num][letter_num][2] = max


class ForwardProfileAligner(Alinger):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        

        
    def process_M(self, state_num, letter_num):
        MM = self.trans_matrix[state_num-1]["MM"] * np.exp(self.res[state_num-1][letter_num-1][0])
        IM = self.trans_matrix[state_num-1]["IM"] * np.exp(self.res[state_num-1][letter_num-1][1])

        if state_num == 1:
            DM = 0
        else:
            DM = self.trans_matrix[state_num-1]["DM"] * np.exp(self.res[state_num-1][letter_num-1][2])
        
        letter = self.seq[letter_num]
        self.res[state_num][letter_num][0] = np.log(self.gen_match[state_num][letter] / self.q[letter]) + np.log(MM + IM + DM)

    def process_I(self, state_num, letter_num):
        if state_num == 0:
            DI = 0
        else:
            DI = self.trans_matrix[state_num]["DI"] * np.exp(self.res[state_num][letter_num-1][2])

        MI = self.trans_matrix[state_num]["MI"] * np.exp(self.res[state_num][letter_num-1][0])
        II = self.trans_matrix[state_num]["II"] * np.exp(self.res[state_num][letter_num-1][1])

        letter = self.seq[letter_num]
        self.res[state_num][letter_num][1] = np.log(self.gen_insert[state_num][letter] / self.q[letter]) + np.log(MI + II + DI)


    def process_D(self, state_num, letter_num):
        if state_num == 1 :
            DD = 0
        else:
            DD = self.trans_matrix[state_num-1]["DD"] * np.exp(self.res[state_num-1][letter_num][2])

        MD = self.trans_matrix[state_num-1]["MD"] * np.exp(self.res[state_num-1][letter_num][0])
        ID = self.trans_matrix[state_num-1]["ID"] * np.exp(self.res[state_num-1][letter_num][1])

        self.res[state_num][letter_num][2] = np.log(MD + ID + DD) 
    
    def fit(self):
        super().fit()
        DM = self.trans_matrix[-1]["DM"] * np.exp(self.res[-1][-1][2])
        MM = self.trans_matrix[-1]["IM"] * np.exp(self.res[-1][-1][1])
        IM = self.trans_matrix[-1]["MM"] * np.exp(self.res[-1][-1][0])
        self.score = DM + MM + IM

