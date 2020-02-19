import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
sys.path.insert(0,currentdir + "/parameter_definitions")

from main_parameter_definitions import *
#from main_parameter_sets_e import *
from main_parameter_sets_ee import *
#from main_parameter_sets_e_fp import *
#from hf_parameter_definitions import *
from challenge_parameter_definitions import *

r5paramSets = [
               ################
               # release d
               ################
               
               r5nd_1kem_0,
               r5nd_3kem_0,
               r5nd_5kem_0,
               r5nd_1pke_0,
               r5nd_3pke_0,
               r5nd_5pke_0,

               r5nd_1kem_5,
               r5nd_3kem_5,
               r5nd_5kem_5,
               r5nd_1pke_5,
               r5nd_3pke_5,
               r5nd_5pke_5,

               r5n1_1kem_0,
               r5n1_3kem_0,
               r5n1_5kem_0,
               r5n1_1pke_0,
               r5n1_3pke_0,
               r5n1_5pke_0,
#
#
#               # Specific use cases
               r5nd_0kem_2iot ,
               r5nd_1kem_4longkey,
               r5n1_3pke_0smallct,

##
##               ################
##               # release e
##               ################
##
               r5nd_1kem_0e,
               r5nd_3kem_0e,
               r5nd_5kem_0e,
#               r5nd_1pke_0e,
#               r5nd_3pke_0e,
#               r5nd_5pke_0e,

               r5nd_1kem_5e,
               r5nd_3kem_5e,
               r5nd_5kem_5e,
#               r5nd_1pke_5e,
#               r5nd_3pke_5e,
#               r5nd_5pke_5e,

               r5n1_1kem_0e,
               r5n1_3kem_0e,
               r5n1_5kem_0e,
#               r5n1_1pke_0e,
#               r5n1_3pke_0e,
#               r5n1_5pke_0e,
               
               # Specific use cases
               r5nd_0kem_2eiot ,
               r5nd_1kem_4elongkey,
               r5n1_3pke_0esmallct,
#
#
#               ################
#               # CHALLENGES
#               ################
#
##               r5n1_1pke_0_challenge_toy,
##               r5n1_1pke_0_challenge_small,
##               r5n1_1pke_0_challenge_medium,
##               #r5n1_1pke_0,
##
##               r5nd_1pke_0_challenge_toy,
##               r5nd_1pke_0_challenge_small,
##               r5nd_1pke_0_challenge_medium,
##               #r5nd_1pke_0,
##
##               r5nd_1pke_5_challenge_toy,
##               r5nd_1pke_5_challenge_small,
##               r5nd_1pke_5_challenge_medium,
##               #r5nd_1pke_5,
##
##               r5nd_1pke_5_challenge_toy_HF,
##               r5nd_1pke_5_challenge_small_HF,
##               r5nd_1pke_5_challenge_medium_HF
##
               
#               ######
#               r5n1_1pke_0e_fp,
#               r5n1_3pke_0e_fp,
               r5n1_5pke_0e_fp,
               r5nd_1pke_0e_fp,
               r5nd_3pke_0e_fp,
               r5nd_5pke_0e_fp,
               r5nd_1pke_5e_fp,
#               r5nd_3pke_5e_fp,
#               r5nd_5pke_5e_fp,
               
               r5nd_3pke_5e_fp4s,
               r5nd_5pke_5e_fp4s,
               r5n1_1pke_0e_fp4f,
               r5n1_3pke_0e_fp4f,
#               r5n1_5pke_0e_fp4f
               

               ####################################
               # SK Distributions paper: 
               # Fixed-weight Ternary parameter set
               ####################################
#               skdistro_fwter_3pke
               
        
               ]

def getParameterSets():
    return r5paramSets
