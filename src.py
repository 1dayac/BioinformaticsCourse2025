from pysuffixarray.core import SuffixArray


class BWTSearcher():
    def __init__(self, text = 'ACAACG'):
        self.text = text; self.n = len(text)
        sufar = SuffixArray(text);      
        self.SA = sufar.suffix_array()

        self.BWT = self.construct_bwt_from_suffix_array()
        self.first_pos = self.get_first_pos_in_SA()
        self.cnt = self.get_cnt_in_BWT()

    def construct_bwt_from_suffix_array(self):
        string, SA = self.text, self.SA

        bwt = []
        n = len(string)
        for i in range(n+1):
            if SA[i] == 0:          # -- the entrance in SA which is the initial string, it's last char = $
                bwt.append('$')
            else:
                bwt.append(string[SA[i] - 1])
        return ''.join(bwt)

    def get_first_pos_in_SA(self):
            string, SA = self.text, self.SA
            first_pos = dict()
            # print(string, len(string))

            for i in range(1, self.n+1): #skipping first $ row
                # print(SA[i])
                char = string[SA[i]]

                if char in first_pos:
                    pass
                else:
                    first_pos[char] = i #first occurence in SA
                    # print(i, string[SA[i]])

            # print(f"text = {string}")
            # print(f"SA_left_column = {SA}")        

            return first_pos

    def get_cnt_in_BWT(self): #, nuc, term):
        # BWT="TC$AAACG"
        BWT = self.BWT
        n = len(BWT)
        
        cnt = {"A":[0]*(n+1),"C":[0]*(n+1),"G":[0]*(n+1),"T":[0]*(n+1) } 
                                            #start with 0 in because we need to discard previous row(first too) in FMindex

        for i in range(1,n+1):
            current_char = BWT[i-1]
                    
            if i > 0: # Copy previous counts
                for key in cnt:
                    cnt[key][i] = cnt[key][i - 1]
            if current_char == "$": #that won't add anything
                continue
            
            cnt[current_char][i] += 1
        
        # print(f"BWT = {list(BWT)}")
        # for k,i in cnt.items():
        #     print(k,i)
        
        return cnt
    

    def check_entrance_in_genome(self, string):
        print(f"searching for '{string}'")
        L, R = 0, self.n 
        # print(self.n), print(len(self.SA))  6 7
        for char in reversed(string):
            if char in self.first_pos:
                first = self.first_pos[char]
            else:
                print(f"Char '{char}' not found in read."); return False
            L = first + self.cnt[char][L] #skipping previous occurences
            R = first + self.cnt[char][R+1] - 1 
                         #R+1 count all occ up to R, -1 to keep [L,R] and not [L,R)
            if L > R: print("No matches found."); return False
            # print(L,R)
        
        idx = []
        for i in range(L,R+1):
            idx.append(str(self.SA[i]))
        print(f"Found read at idx: {','.join((reversed(idx)))}.")
        # return reversed(idx)
        return True


# textS = ['ACAACG', 'ACGT', 'AAAAAAAAAAAAAAAAAA']
# # textS = ['ACAACG']
# for k, text in enumerate(textS):
#     print(f"{k}-th test, text = {text}, last index = {len(text)-1}")
#     searcher = BWTSearcher(text=text)
#     # print(searcher.BWT)
#     # print(len(searcher.BWT)==len(searcher.SA), len(searcher.text)==len(searcher.SA) - 1) # 7 7 6
#     # print(searcher.first_pos)
#     for pattern in ["CA", "ACA", "T", "AA"]:
#         searcher.find_idx_of_string_in_text(pattern)

#     print(f"SUCCESS on {k}st text {text}\n")
        
# print(searcher.find_idx_of_string_in_text('CA'))








