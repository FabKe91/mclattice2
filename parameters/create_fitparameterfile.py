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
PLCPAIRS = ["DPPC_CHL1", "DUPC_CHL1", "DUPC_CHL1"]

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
    print(df)
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
    print(df)
    growth = df.growth.iloc[0]
    deltaN = df.deltaN.iloc[0]
    Nmax = df.Nmax.iloc[0]
    transitionT = df.transitionT
    return ' '.join([str(Emax), str(Egain), str(growth), str(S0) ])

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
    print(df)
    deg = df.deg.iloc[0] #=k
    E0 = df.E0.iloc[0]
    S0 = df.S0.iloc[0]
    Sshift = 1.5
    return ' '.join([str(S0), str(Sshift), str(deg), str(E0) ])

def _read_poly_parafile(df):
    ''' for polynomials of the for \sumi=0 coeff[i]*x**i '''
    print("Warning: polynomial coefficients must be given in ascending order")
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

    sep = "\t"

    try:
        dat = pd.read_csv(parafilename, sep=sep)
    except FileNotFoundError(f"{parafilename} not found"):
        print("Missing parameter file for Eplpl(S, Nc) interaction.")
        raise
    
    if len(dat.columns) == 1:
        print("Failed to read Eplpl file")
        raise ValueError(f"Separator {sep} incorrect")

    print(dat)
    for plpair in PLPLPAIRS:
        for Nc in range(7):
            df = dat.loc[(dat.cholconc == CHOLCONC) & (dat.pair == plpair) & (dat.Nc == Nc)]
            print(df)
            outstr = ' '.join(["Enthalpy", plpair.replace("_"," "), str(Nc), _read_logistic_parafile(df), "\n"]) 
            f.write(outstr)


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




def add_Ecc(f, parafilename="fitparameters_CC_Nc.txt"):
    ''' linear fit '''
    sep = "\t"

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


def add_selfE(f, parafilename="fitparameters_selfE.txt"):
    ''' polynomial fit '''
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


def add_NLL():
    ''' '''

def add_NLC():
    ''' '''

def add_entropy():
    ''' '''


###############################



if __name__ == "__main__":

    with open(OUTFNAME, "w") as f:
        f.write(f"# Parameters from c{CHOLCONC} and T{TEMPERATURE}\n")

        add_Eplpl(f)
        add_Eplc(f)
        #add_Ecc(f)
        #add_selfE(f)
        #add_NLL(f)
        #add_NLC(f)
        #add_entropy(f)

    





