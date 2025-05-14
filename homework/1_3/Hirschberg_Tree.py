from anytree import Node as ATNode, RenderTree  # ADD THIS AT THE TOP

class Hirsberg_Tree():
    def __init__(self, Astr, Bstr, ins=-2, dele=-2, match=2, wrong=-1):
        self.n = len(Astr)
        self.m = len(Bstr)
        self.Astr = Astr
        self.Bstr = Bstr
        self.ins = ins
        self.dele = dele
        self.match = match
        self.wrong = wrong
        self.tree = self.Hirsberg_align(self.Astr, self.Bstr)

    class node():
        def __init__(self, LEFTstr=None, Rstr=None):
            self.left, self.right = LEFTstr, Rstr
            self.left_child = None
            self.right_child = None

    def scoring_row(self, X, Y):       #Linear time , only keep 2 rows while counting
        gap = -1
        match = 1
        wrong = -1
        n, m = len(X), len(Y)
        old_mat = [0]*(m+1)
        new_mat = [0]*(m+1)

        for j in range(m + 1):
            old_mat[j] = self.ins * j
        
        if n == 0:
             return old_mat
        
        for i in range(1, n + 1):
            new_mat = [0]*(m+1)
            new_mat[0] = self.dele * i #same init, 1st col differs

            for j in range(1, m + 1):
                if X[i - 1] == Y[j - 1]:  #left                          #top_prev              #top-diag
                    new_mat[j] = max(new_mat[j - 1] + self.ins, old_mat[j] + self.dele, old_mat[j - 1] + self.match)
                else: 
                    new_mat[j] = max(new_mat[j - 1] + self.ins, old_mat[j] + self.dele, old_mat[j - 1] + self.wrong)
            old_mat = new_mat.copy()
            # print(old_mat)
            # print(new_mat)

        return new_mat ## [n] only last row needed
    

    def Hirsberg_align(self, X=None, Y=None):
        if X is None and Y is None:
            X, Y = self.Astr, self.Bstr

        if X == "":
            return self.node("-" * len(Y), Y)
        if Y == "":
            return self.node(X, "-" * len(X))

        tree = self.node(X, Y)

        if len(X) == 1 or len(Y) == 1: #leafy leaf
            return self.node(X, Y)

        midX = len(X) // 2
        score1 = self.scoring_row(X[:midX], Y)
        score2 = self.scoring_row(X[midX:][::-1], Y[::-1]) #this way they meet in the middle
        score2.reverse()

        for i in range(len(Y)+1):
            score1[i] += score2[i]
        y_split = argmax(score1)

        tree.left_child = self.Hirsberg_align(X[:midX], Y[:y_split])
        tree.right_child = self.Hirsberg_align(X[midX:], Y[y_split:])
        return tree

    def get_aligned_sequences(self, tree):
        if tree is None:
            return ("", "")
        
        if tree.left_child is None and tree.right_child is None:
            return (tree.left, tree.right)
        
        left_A, left_B = self.get_aligned_sequences(tree.left_child)
        right_A, right_B = self.get_aligned_sequences(tree.right_child)

        return (left_A + right_A, left_B + right_B)

    # ======== NEW METHOD HERE =========
    def tree_terminal(self):        #Visualisation took from GPT
        def to_anytree(n, parent=None):
            label = f"({n.left},{n.right})"
            at = ATNode(label, parent=parent)
            if n.left_child:
                to_anytree(n.left_child, at)
            if n.right_child:
                to_anytree(n.right_child, at)
            return at
        
        at_root = to_anytree(self.tree)
        for pre, _, nd in RenderTree(at_root):
            print(f"{pre}{nd.name}")
    # ===================================


def argmax(lis):
    m = - float("inf")
    ans = 0
    for i, x in enumerate(lis):
        if x > m:
            m = x
            ans = i
    return ans


def main():
    htree = Hirsberg_Tree("AGTACGCA", "TATGC")
    
    # htree = Hirsberg_Tree("xyzAGTACGCA", "AGCTGASasdas")

    htree.tree_terminal()

    aligned_A, aligned_B = htree.get_aligned_sequences(htree.tree)
    print("\nВыравнивание последовательностей (Alignment):")
    print(aligned_A)
    print(aligned_B)

main()
