import sys
import pandas as pd 

if len(sys.argv) != 3+1:
    print("Input arguments required: <CHOLCONC> <T> <outputfilename>")
    sys.exit()

CHOLCONC    = int(sys.argv[1])
TEMPERATURE = int(sys.argv[2])
OUTFNAME = sys.argv[3]

PLLIST = ["DPPC", "DUPC"]
CHOLNAME = "CHL1"

PLPLPAIRS = ["DPPC_DPPC", "DUPC_DUPC", "DUPC_DPPC"]
PLCPAIRS = ["DPPC_CHL1", "DUPC_CHL1"]

#######################################

def _read_logistic_parafile(df):
    '''
     /**
     * NOTE: Primarily used for PL-PL interactions (function of order)
     * NOTE: This function is similar to sigmoid, but introduced to discern from sigmoid used for neighbor scaling function
     * logistic function of form:
     *    f(S) = Emax + ( Egain / (1+exp(-k(S-S0))) )
     * @param[in] coeff list of coeffs: <Emax> <Egain> <k> <S0>
     * @param[in] S function variable
     */
    '''
    growth = df.growth.iloc[0]
    Egain = df.Egain.iloc[0]
    Emax = df.Emax.iloc[0]
    S0 = 1.0
    return ' '.join([str(Emax), str(Egain), str(growth), str(S0) ])

def _read_sigmoid_parafile(df):
    '''
        /**
         * NOTE: Primarily used for neighbor rescaling of DPPC (function of temperature)
         * NOTE: This function is similar to sigmoid, but introduced to discern from sigmoid used for neighbor scaling function
         * logistic function of form:
         *    f(S) = Nmax + ( deltaN / (1+exp(-k(T-transitionT))) )
         * @param[in] coeff list of coeffs: <Nmax> <deltaN> <k> <transitionT>
         * @param[in] T function variable
         */
    '''
    growth = df.growth.iloc[0]
    deltaN = df.deltaN.iloc[0]
    Nmax = df.Nmax.iloc[0]
    transitionT = df.transitionT.iloc[0]
    return ' '.join([str(Nmax), str(deltaN), str(growth), str(transitionT) ])

def _read_lennard_parafile(df):
    '''
     /**
     * NOTE: Primarily used for PL-CHOL interaction
     * lennard jones function of form:
     *      Q = [(S0 - Sshift) / (S - Sshift)]**k
     *      f(S) = -E0 * Q * (Q - 2)
     * @param[in] coeff list of coeffs: <S0> <Sshift> <k> <E0>
     * @param[in] S function variable
     */
    '''
    deg = df.deg.iloc[0] #=k
    E0 = df.E0.iloc[0]
    S0 = df.S0.iloc[0]
    Sshift = 1.5
    return ' '.join([str(S0), str(Sshift), str(deg), str(E0) ])

def _read_poly_parafile(df):
    ''' for polynomials of the for \sumi=0 coeff[i]*x**i '''
    print("Warning: polynomial coefficients must be given in ascending order")
    print("data okay?\n", df)
    return ' '.join([str(df[col].iloc[0]) for col in df.columns])




### ##########################

def add_Eplpl(f, parafilename="fitparameters_EPLPL_Nc.txt"):
    '''
     /**
     * NOTE: Primarily used for PL-PL interactions (function of order)
     * NOTE: This function is similar to sigmoid, but introduced to discern from sigmoid used for neighbor scaling function
     * logistic function of form:
     *    f(S) = Emax + ( Egain / (1+exp(-k(S-S0))) )
     * @param[in] coeff list of coeffs: <Emax> <Egain> <k> <S0>
     * @param[in] S function variable
     */
    '''
    print("reading: ", parafilename)

    f.write("#logistic function - Eplpl(Smin, Nc):\n")
    f.write("# f(S) = Emax + ( Egain / (1+exp(-k(S-S0))) )\n")
    f.write("# <Emax> <Egain> <k> <S0>\n")


    sep = "\t"

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for Eplpl(S, Nc) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read Eplpl file")
        raise ValueError(f"Separator {sep} incorrect")

    for plpair in PLPLPAIRS:
        for Nc in range(7):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.pair == plpair) & (dat.Nc == Nc)]
            outstr = ' '.join(["Enthalpy", plpair.replace("_"," "), str(Nc), _read_logistic_parafile(df), "\n"]) 
            f.write(outstr)

    f.write("\n")

def add_Eplc(f, parafilename="fitparameters_EPLC_Nc.txt"):
    '''
     /**
     * NOTE: Primarily used for PL-CHOL interaction
     * lennard jones function of form:
     *      Q = [(S0 - Sshift) / (S - Sshift)]**k
     *      f(S) = -E0 * Q * (Q - 2)
     * @param[in] coeff list of coeffs: <S0> <Sshift> <k> <E0>
     * @param[in] S function variable
     */
    '''
    sep = "\t"

    f.write("#lennard jones like - Eplc(Shost, Nc):\n")
    f.write("# Q = [(S0 - Sshift) / (S - Sshift)]**k\n")
    f.write("# f(S) = -E0 * Q * (Q - 2)\n")
    f.write("#\t\t\t<S0> <Sshift> <k> <E0>\n")

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for Eplc(Spl, Nc) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read Eplc file")
        raise ValueError(f"Separator {sep} incorrect")

    for plcpair in PLCPAIRS:
        for Nc in range(6):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.pair == plcpair) & (dat.Nc == Nc)]
            outstr = ' '.join(["LipidCholEnergie", plcpair.replace(f"_{CHOLNAME}", ""), str(Nc), _read_lennard_parafile(df), "\n"]) 
            f.write(outstr)
    f.write("\n")




def add_Ecc(f, parafilename="fitparameters_ECC_Nc.txt"):
    ''' linear fit '''
    sep = "\t"
    print("reading: ", parafilename)

    f.write("#linear - Ecc(Nc):\n")

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for Ecc(Nc) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read Ecc file")
        raise ValueError(f"Separator {sep} incorrect")

    for pl in PLLIST:
        df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.PL == pl)].drop(columns=["cholconc", "PL"])
        outstr = ' '.join(["CholCholEnergie", pl, _read_poly_parafile(df), "\n"]) 
        f.write(outstr)
    f.write("\n")


def add_selfE(f, parafilename="fitparameters_selfE.txt"):
    ''' polynomial fit '''
    print("reading: ", parafilename)

    f.write("#polynomial - self interaction Eself(S):\n")

    sep = "\t"

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for Eself(S) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read selfE file")
        raise ValueError(f"Separator {sep} incorrect")

    for pl in PLLIST:
        df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.PL == pl)].drop(columns=["cholconc", "PL"])
        outstr = ' '.join(["SelfEnergie", pl, _read_poly_parafile(df), "\n"]) 
        f.write(outstr)

    f.write("\n")

def add_NLL(f, lipidtype, parafilename="fitparameters_NofT_Nc_{}.txt"):
    ''' fit depends on lipid type '''

    parafilename = parafilename.format(lipidtype)
    print("reading: ", parafilename)

    f.write("#Lipid-Lipid neighbor functions sigmoidal(DPPC)/polynom(DUPC) N_pl(T, Nc):\n")
    if lipidtype == "DPPC":
        f.write("# f(S) = Nmax + ( deltaN / (1+exp(-k(T-transitionT))) )\n")
        f.write("#\t\t\t <Nmax> <deltaN> <k> <transitionT>\n")

    sep = "\t"

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for NLL(T, Nc) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read NLL(T, Nc) file")
        raise ValueError(f"Separator {sep} incorrect")

    if lipidtype == "DPPC":
        for Nc in range(5):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.PL == lipidtype) & (dat.Nc == Nc)].drop(columns=["cholconc", "PL", "Nc"])
            outstr = ' '.join(["LipidLipidNeighPara", lipidtype, _read_sigmoid_parafile(df), "\n"]) 
            f.write(outstr)
    elif lipidtype == "DUPC":
        for Nc in range(5):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.PL == lipidtype) & (dat.Nc == Nc)].drop(columns=["cholconc", "PL", "Nc"])
            outstr = ' '.join(["LipidLipidNeighPara", lipidtype, _read_poly_parafile(df), "\n"]) 
            f.write(outstr)
    else:
        raise ValueError(f"Unknown lipidtype: {lipidtype}")
    f.write("\n")


def add_NLC(f, lipidtype, parafilename="fitparameters_NofT_Nc_CHOL-{}.txt"):
 
    parafilename = parafilename.format(lipidtype)
    print("reading: ", parafilename)

    f.write("#Chol-Lipid neighbor functions sigmoidal(DPPC)/polynom(DUPC) N_pl(T, Nc):\n")
    if lipidtype == "DPPC":
        f.write("# f(S) = Nmax + ( deltaN / (1+exp(-k(T-transitionT))) )\n")
        f.write("#\t\t\t <Nmax> <deltaN> <k> <transitionT>\n")


    sep = "\t"

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for NLL(T, Nc) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read NLL(T, Nc) file")
        raise ValueError(f"Separator {sep} incorrect")

    if lipidtype == "DPPC":
        for Nc in range(5):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.PL == lipidtype) & (dat.Nc == Nc)].drop(columns=["cholconc", "PL", "Nc"])
            outstr = ' '.join(["CholLipidNeigh", lipidtype, _read_sigmoid_parafile(df), "\n"]) 
            f.write(outstr)
    elif lipidtype == "DUPC":
        for Nc in range(5):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.PL == lipidtype) & (dat.Nc == Nc)].drop(columns=["cholconc", "PL", "Nc"])
            outstr = ' '.join(["CholLipidNeigh", lipidtype, _read_poly_parafile(df), "\n"]) 
            f.write(outstr)
    else:
        raise ValueError(f"Unknown lipidtype: {lipidtype}")
    f.write("\n")


def add_entropy(f, parafilename="fitparameters_entropies.txt"):
    ''' '''

def add_pscd(f, parafilename="fitparameter_pscd.txt"):
    ''' '''


###############################



if __name__ == "__main__":

    with open(OUTFNAME, "w") as f:
        f.write(f"# Parameters from c{CHOLCONC} and T{TEMPERATURE}\n\n")

        add_Eplpl(f)
        add_Eplc(f)
        add_Ecc(f)
        add_selfE(f)
        add_NLL(f, "DPPC")
        add_NLC(f, "DPPC")
        add_NLL(f, "DUPC")
        add_NLC(f, "DUPC")
        #add_entropy(f)
        #add_pscd(f)

        f.write("\n")

    





