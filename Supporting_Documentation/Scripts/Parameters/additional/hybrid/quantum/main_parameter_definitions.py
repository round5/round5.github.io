# generated file

import sys, inspect, os
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
from r5_parameter_set import R5_paramSet

# Round5 ring parameter sets without error correction.
r5nd_1kem_0 = R5_paramSet("R5ND_1KEM_0d",618,618,104,11,8,4,1,1,1,16,0,0,False)
r5nd_3kem_0 = R5_paramSet("R5ND_3KEM_0d",786,786,384,13,9,4,1,1,1,24,0,0,False)
r5nd_5kem_0 = R5_paramSet("R5ND_5KEM_0d",1018,1018,428,14,9,4,1,1,1,32,0,0,False)

r5nd_1pke_0 = R5_paramSet("R5ND_1PKE_0d",586,586,182,13,9,4,1,1,1,16,0,0,False)
r5nd_3pke_0 = R5_paramSet("R5ND_3PKE_0d",852,852,212,12,9,5,1,1,1,24,0,0,False)
r5nd_5pke_0 = R5_paramSet("R5ND_5PKE_0d",1170,1170,222,13,9,5,1,1,1,32,0,0,False)

# Round5 ring parameter sets with error correction.
r5nd_1kem_5 = R5_paramSet("R5ND_1KEM_5d",490,490,162,10,7,3,1,1,1,16,5,190,True)
r5nd_3kem_5 = R5_paramSet("R5ND_3KEM_5d",756,756,242,12,8,2,1,1,1,24,5,218,True)
r5nd_5kem_5 = R5_paramSet("R5ND_5KEM_5d",940,940,414,12,8,2,1,1,1,32,5,234,True)

r5nd_1pke_5 = R5_paramSet("R5ND_1PKE_5d",508,508,136,10,7,4,1,1,1,16,5,190,True)
r5nd_3pke_5 = R5_paramSet("R5ND_3PKE_5d",756,756,242,12,8,3,1,1,1,24,5,218,True)
r5nd_5pke_5 = R5_paramSet("R5ND_5PKE_5d",946,946,388,11,8,5,1,1,1,32,5,234,True)

# Round5 non-ring parameter sets
# Note: setting the "switchrings" parameter=True for non-ring variants
# due to absence of correlation issue associated with reductions
# modulo the prime-order cyclotomic polynomial
r5n1_1kem_0 = R5_paramSet("R5N1_1KEM_0d",594,1,238,13,10,7,3,7,7,16,0,0,True)
r5n1_3kem_0 = R5_paramSet("R5N1_3KEM_0d",881,1,238,13,10,7,3,8,8,24,0,0,True)
r5n1_5kem_0 = R5_paramSet("R5N1_5KEM_0d",1186,1,712,15,12,7,4,8,8,32,0,0,True)

r5n1_1pke_0 = R5_paramSet("R5N1_1PKE_0d",636,1,114,12,9,6,2,8,8,16,0,0,True)
r5n1_3pke_0 = R5_paramSet("R5N1_3PKE_0d",876,1,446,15,11,7,3,8,8,24,0,0,True)
r5n1_5pke_0 = R5_paramSet("R5N1_5PKE_0d",1217,1,462,15,12,9,4,8,8,32,0,0,True)

# Round5 parameters for specific use cases
r5nd_0kem_2iot     = R5_paramSet("R5ND_0KEM_2iot",372,372,178,11,7,3,1,1,1,16,2,53,True)
r5nd_1kem_4longkey = R5_paramSet("R5ND_1KEM_4longkey",490,490,162,10,7,3,1,1,1,24,4,163,True)
r5n1_3pke_0smallct = R5_paramSet("R5N1_3PKE_0smallCT",757,1,378,14,9,4,1,192,1,24,0,0,True)

# Some other alternatives that are not submitted, but that could be interesting.
r5nd_5kem_5td   = R5_paramSet("R5ND_5KEM_5td",946,946,388,12,8,2,1,1,1,32,5,234,True)
r5nd_1pke_5td   = R5_paramSet("R5ND_1PKE_5td",522,522,124,10,7,4,1,1,1,16,5,190,True)
r5nd_5pke_5td   = R5_paramSet("R5ND_5PKE_5td",946,946,388,12,8,3,1,1,1,32,5,234,True)
r5nd_0kem_3iots = R5_paramSet("R5ND_0KEM_3iots",372,372,178,11,7,3,1,1,1,16,3,91,True)
r5nd_0pke_4iot  = R5_paramSet("R5ND_0PKE_4iot",418,418,92,11,7,3,1,1,1,16,4,149,True)
