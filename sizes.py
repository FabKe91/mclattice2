import seaborn as sns 

def set_context():
    sns.set_context("poster")
    #sns.set_context("paper", font_scale=1.2, rc={"line.linewidth":1, **rc})


cm = 1/2.54

a169 = 1/(16/9) 
a1610 = 1/(16/10)
a43  =  1/(4/3)
a54  =  1/(5/4)

singlex=16*cm

rc = {
        "ytick.major.size"     : 5*0.6 ,   ## major tick size in points
        "ytick.minor.size"     : 3.5*0.6   ,   ## minor tick size in points
        "ytick.major.width"    : 1.3*0.6 ,   ## major tick width in points
        "ytick.minor.width"    : 0.9*0.6 ,   ## minor tick width in points
        "xtick.major.size"     : 5*0.6 ,   ## major tick size in points
        "xtick.minor.size"     : 3.5*0.6   ,   ## minor tick size in points
        "xtick.major.width"    : 1.3*0.6 ,   ## major tick width in points
        "xtick.minor.width"    : 0.9*0.6 ,   ## minor tick width in points
        }






