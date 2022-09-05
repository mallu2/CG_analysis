#!/Users/mallu899/anaconda3/bin/python


class Iterator_Nb:
    """Iterate over coordinates--> adjust for each model: monomer: 59
                                                    dimer: 118
                                                    N-core: 288
                                                    full-lenth: 
    Define the number of models and set step.
    Retruns the coordinates fo the beads in the trajectory file.
    lp = 288 length of the Ncore sequence,
    lp = 59 length of the monomer sequence"""
    
    def __init__(self):
        self = self
 
    def readlinesfromfile(self):
        """Returns list of lines in the file."""
        lines = []
        with open(self, "r") as f:
            lines.append(f.readlines())
        return(lines)
    
    def traj(file,framecount,step):
        lines = Iterator_Nb.readlinesfromfile(file)
        lp = 288
        mods = []
        for i in range(0,framecount,step):
            mods.append(lines[0][(lp+607)+i*(lp+2):(lp+607)+(i+1)*(lp+2)][1:(lp+1)])
        return(mods)
    

    def pdb(file):
        """Retruns coordinates of CA atoms in pdb."""
        lines = Iterator_N.readlinesfromfile(file)
        mods = []
        for i in range(len(lines)):
            if 'CA' in lines[i]:
                mods.append(lines[i])
        return(mods)
class Region():
    """Retruns coordinates of distinct regions in lacI protein."""

    def H1A(data): 
        H1A = []
        for frame in data:
            H1A.append(frame[5:13])
        return(H1A)

    def H1B(data):
        lp = len(data[0])
        mono=lp/2 
        H1B = []
        B1 = int(5+mono)
        B2 = int(13+mono)
        for frame in data:
            H1B.append(frame[B1:B2])
        return(H1B)
    
    def H2A(data): 
        H2A = []
        for frame in data:
            H2A.append(frame[16:25])
        return(H2A)

    def H2B(data):
        lp = len(data[0])
        mono=lp/2 
        H2B = []
        B1 = int(16+mono)
        B2 = int(25+mono)
        for frame in data:
            H2B.append(frame[B1:B2])
        return(H2B)
    
    def H3A(data): 
        H3A = []
        for frame in data:
            H3A.append(frame[31:43])
        return(H3A)

    def H3B(data):
        lp = len(data[0])
        mono=lp/2 
        H3B = []
        B1 = int(31+mono)
        B2 = int(43+mono)
        for frame in data:
            H3B.append(frame[B1:B2])
        return(H3B)

    def hiA(data):
        hiA = []
        for frame in data:
            hiA.append(frame[50:56])
        return(hiA)

    def hiB(data):
        lp = len(data[0])
        mono=lp/2 
        hiB = []
        B1 = int(50+mono)
        B2 = int(56+mono)
        for frame in data:
            hiB.append(frame[B1:B2])
        return(hiB)
