# common parameters

from itertools import count
KEYFHOF = "FINAL HEAT OF FORMATION"
KEYCA = "COSMO AREA"
KEYDEFAULT = "PUT KEYWORDS HERE"
KEYCOSMO = "NR. ELEM.             COORDINATES     RADIUS  COSMO-CHARGE" #to position cursor in COSMO file
MOPACPATH = "/opt/mopac/MOPAC2012_64.exe"
DEFCORE = 4
DEFPARAMFILE = "list_param.in"
DEFREFFILE = "list_dgexp.in"
DEFREPORFILE = "report.log"
DEFSUMMARYFILE = "summaryfile.log"
DEFTEMPLATEDIR= "templates"
DEFINPUTEXT = ".mop"
DEFTEMPLATEEXT = ".mop"
DEFGASDIR = "gas"
DEFOPTIMIZEPARAM = "optimize.log"
RCONST=8.31 # gases constant in J K-1 mol-1
TEMP = 298 # temperature in K
NS = 3.33679E+028 # numeral density in m-3
PI = 3.14159
J2KCAL = 0.0002388458966275
step = ("%04i" % i for i in count(1)) #work with itertools to make 002, 003 format name
CONFIG_FILENAME = 'config.in'

