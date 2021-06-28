############################################################
#STANDARD CENTURY (MATRIX)
#It simulates the C dynamics over the experiment length
############################################################
# matrix representation of Century/ORCHIDEE, 7 pools;
# aboveground metabolic litter; belowground meta litter; above structure
# litter; below structure litter; active SOC; slow SOC; passive SOC
# Yuanyuan Huang yuanyuanhuang2011@gmail.com
############################################################
#translation to Python
#Elisa Bruni ebruni93@gmail.com
############################################################
import sys
import numpy as npy
import scipy
from scipy.optimize import minimize
from scipy.optimize import least_squares
import numdifftools as ndt
import math
from random import gauss
import xlrd
import pandas as pd
import time
import datetime
from datetime import datetime

npy.set_printoptions(threshold=sys.maxsize)
#npy.set_printoptions(linewidth=npy.inf)


out=open("SOC_data_model_10.txt","wb")

##################
#open optimized struc:metab ratios
#################
n=0 # contatore
inp=open('opt_abfractions_forscript2.9.txt','rb')
while inp.read(1):
    inp.seek(-1,1)
    XXmet=npy.load(inp)
n+=1


ns=XXmet.shape[0]
XXstr=1-XXmet

imp_frac = npy.zeros((ns,4))
imp_frac[:,0]=XXstr[:,0]
imp_frac[:,1]=XXmet[:,0]
imp_frac[:,2]=XXstr[:,1]
imp_frac[:,3]=XXmet[:,1]



#########################

def water_rh(clay,water_in_m3):

        ######################################
        #conversion of water_in to relative humidity [coeff (0;1)]
        #####################################
        # water_in= (mcs-mcw)/(mcfc-mcw)
        # mcs = orchidee outputs (soil hum = mrsos)
        # mcw = welting point
        # mcfc = field capacity

        #default values for orchidee (coarse, medium, fine)
        coarse=0
        medium=1
        fine=2

        mcfc=npy.array([0.1218, 0.1654, 0.2697])
        mcw= npy.array([0.0657,  0.0884, 0.1496])

        #relative humidity (water_in) is a value between 0 and 1
        # if water_in >1 -> it means water content is above field capacity
        # -> take the minimum between 1 and the actual value

        option_1=False

        ##########################################
        #OPTION 1
        ##choose one type of soil separate between coarse, medium and fine (sandy->coarse; silty->medium; clay->fine)
        ##########################################
        if(option_1):
                site_texture=coarse
                print"SOIL TYPE:",site_texture

                water_in = npy.minimum(1.0,(map(float,water_in_m3) - mcw[site_texture])/(mcfc[site_texture]-mcw[site_texture]))
	        water_in = npy.maximum(0.0,water_in)
                print "mean water_in",npy.mean(water_in)

        ##########################################
        #OPTION 2
        #Weighted water_in
        ##########################################
        else:
                silt = (1-clay)/2 #suppose %silt=%sand ->see French soil map
                sandy= (1-clay)/2

                weighted_water_in = clay*(map(float,water_in_m3) - mcw[fine])/(mcfc[fine]-mcw[fine])+silt*((map(float,water_in_m3) - mcw[medium])/(mcfc[medium]-mcw[medium]))+sandy*((map(float,water_in_m3) - mcw[coarse])/(mcfc[coarse]-mcw[coarse]))
                water_in = npy.minimum(1.0,weighted_water_in)
		water_in = npy.maximum(0.0,water_in)
                print "mean water_in_m3",npy.mean(water_in_m3)
                print "mean water in",npy.mean(water_in)
	return water_in

#########################################################
def AB_NanZeroRemover(site_T0,site_T0_name,iout,ROOTDIR):
#########################################################
        if ( iout > 0 ):
                out1=open(ROOTDIR+"AB.data","wb")
        #######################################################################
        # remove NaNs and ZEROs from ABOVE and BELOW
        #######################################################################
        yy=npy.asarray(site_T0['Year']).astype(npy.int16)             # array of years
        aa=npy.asarray(site_T0['ABOVE']).astype(npy.float64)*100/365  # array of C_above (gC/m2)
        bb=npy.asarray(site_T0['BELOW']).astype(npy.float64)*100/365  # array of C_below (gC/m2)
        aa0=npy.where(npy.isnan(aa),0,aa)  # replace NaN with zero in "above" 
        YEAR=yy[aa0>0]                     # select years where above>0
        abo=aa[aa0>0]                      # select ABOVE>0
        bel=bb[aa0>0]                      # select corresponding BELOW
        if (iout > 0):
                XX=npy.stack((YEAR,abo,bel),axis=0)
                npy.save(out1,XX)
        print site_T0_name,': AB_NanZeroRemover --> selected ',len(YEAR),' out of ',len(yy)
        return abo,bel,YEAR

#############################################
def SOC_NanZeroRemover(site_T0,site_T0_name):
#############################################
        BigRelativeError=0.10
        # put year, soc, variance into numpy arrays
        yy0=npy.asarray(site_T0['Year']).astype(npy.int16)            # array of years
        ss0=npy.asarray(site_T0['SOC']).astype(npy.float64)           # array of SOCs
        vv0=npy.asarray(site_T0['SOC variance']).astype(npy.float64)  # array of SOC variances
        ss0=npy.where(npy.isnan(ss0),0,ss0)      # replace NaN with 0
        sc=ss0[ss0>0]                            # cut away all 0s, sc now corresponds to real measurements
        YEAR=yy0[ss0>0]                          # select the years corresponding to sc
        sc=sc*100                                # pass to gC/m2
        vv0=vv0*10000
        std2=npy.std(sc)**2                      # square standard deviation of the measurements (use when no error provided - ??)
        if (std2 == 0):
                std2=(BigRelativeError*sc)**2    # <-- check  # if std == 0 use BigRelativeError in spite of std
        vv0=npy.where(npy.isnan(vv0),std2,vv0)   # Replace NaN in variance array with std2
        vv0=npy.where(vv0==0,std2,vv0)           # Replace 0 in variance with std2
        var=vv0[ss0>0]                           # Restrict variance corresponding to the selected SOCs data
        print site_T0_name,': SOC_NanZeroRemover (cleanup of SOC data) --> selected ',len(YEAR), ' years out of ',len(yy0)
        return sc,var,YEAR

#################################################
#
# INITIALIZATION
#
#################################################
NEW_ITER = 0
np = 7
one_year = 365
one_day = 86400
dt=1 # daily time step

iforce_recycle=30*one_year

prior_soilQ10 = npy.log(2)
prior_t = 30.
prior_soilQ10_t = npy.array([prior_soilQ10,prior_t])

Q10 = 10.

frac_soil_metab_aa  = 0.45  # aboveground metabolic to active SOC
frac_soil_metab_ab  = 0.45  # below metabolic to active SOC 
frac_soil_struct_aa = 0.55 # above structure to active SOC 
frac_soil_struct_ab = 0.45 # below structure to active SOC
frac_soil_struct_sa = 0.7  # above structure to slow SOC
frac_soil_struct_sb = 0.7  # below structure to slow SOC

frac_passive_active = 0.004 # active to passive
frac_active_slow    = 0.42     # slow to active
frac_passive_slow   = 0.03    # slow to passive
frac_active_passive = 0.45  # passive to active
frac_slow_passive   = 0.0     # passive to slow

lignin_struc_cmatrix = npy.array([0.76, 0.72]) # aboveground lignin in struc litter; belowground lignin in structure litter

tau_metabolic =   0.066*one_year  # turnover time per day
tau_struct    =   0.245*one_year
tau_active    =   0.149*one_year
tau_slow      =   5.480*one_year
tau_passive   =     241*one_year      # try with higher passive tau

prior_tau=npy.array([tau_metabolic,tau_struct,tau_active,tau_slow,tau_passive])

flux_tot_coeff = [1.2, 1.4, 0.75] #only the third is used
litter_struct_coef = 3. 

CHI2_PRINT_FREQUENCY=50


######################################################
#For each site: Set SITE name and experiment duration
#####################################################
#entire dataset
#ROOTDIR=Root_directory
#loc_exp = ROOTDIR+experiment_location

C_input_exp = pd.read_excel(loc_exp)
site_names = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names = map(str, site_names)
print "SITE NAMES", site_names

#Control plots names
site_T0_array=npy.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 'LAJA3_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 'DwN0', 'TREV1_Min','AVRI_T12TR'])

#Stationary solution array for each experiment T0
N_sites=len(site_names)

SOC_exp_array=[]
SOC_clean_exp_array=[]
SOC_clean_exp_variance=[]
SOC_clean_year=[]

SITE_year0=npy.zeros(N_sites)
SITE_date_init=npy.zeros(N_sites)
SITE_date_end=npy.zeros(N_sites)
SITE_date_init_ss=npy.zeros(N_sites)
SITE_date_end_ss=npy.zeros(N_sites)
SITE_exper_len = npy.zeros(N_sites)
SITE_clay=npy.zeros(N_sites)
SITE_water_t=npy.zeros(N_sites)
SITE_temp_t=npy.zeros(N_sites)
SITE_water_in=npy.zeros(N_sites)
SITE_temp_in=npy.zeros(N_sites)
SITE_ABOVE_mean=npy.zeros(N_sites)
SITE_BELOW_mean=npy.zeros(N_sites)
SITE_ERR2_ABOVE_mean=npy.zeros(N_sites)
SITE_ERR2_BELOW_mean=npy.zeros(N_sites)
SITE_mean_relhum=npy.zeros(N_sites)
SITE_litterinc = npy.zeros((N_sites,4))
SITE_litterinc_err2 = npy.zeros((N_sites,4))
SITE_water_in=[]
SITE_temp_in=[]
SITE_water_t=[]
SITE_temp_t=[]
SITE_ABOVE=[]
SITE_BELOW=[]

#Sites
CHNO3=0
COL=1
CREC3=2
FEU=3
JEU1=4
LAJA2=5
LAJA3=6
RHEU1=7
RHEU2=8
ARAZ=9
ULTU=10
BROAD=11
FOGGIA=12
TREV1=13
AVRI=14

# use all sites
Use_Site=npy.arange(N_sites)

#select sites to be used
Use_Site=npy.zeros(N_sites,dtype=npy.int16)
Use_Site[CHNO3]=1
Use_Site[COL]=1
Use_Site[CREC3]=1
Use_Site[FEU]=1
Use_Site[JEU1]=0
Use_Site[LAJA2]=1
Use_Site[LAJA3]=1
Use_Site[RHEU1]=1
Use_Site[RHEU2]=1
Use_Site[ARAZ]=1
Use_Site[ULTU]=1
Use_Site[BROAD]=1
Use_Site[FOGGIA]=1
Use_Site[TREV1]=1
Use_Site[AVRI]=1

j=0
for site in site_names:
        ##########################################
        #IMPORT C input
        ##########################################
        #import metabolic:structural fraction at site
	frac_array=imp_frac[j]

        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        print "READING DATA OF SITE: ",site
        print ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
        site_df = C_input_exp[(C_input_exp['ID.Site'].values == [site])]
        year_0 = npy.min(site_df['Year'])
        year_end = npy.max(site_df['Year'])
        site_T0_name = site_T0_array[j]
        site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]
        SOC_dyn= site_T0['SOC']
        SOC_dyn=SOC_dyn*100 #(gC/m2)
        SOC_dyn = SOC_dyn.astype(float).interpolate()
        SOC_dyn = pd.DataFrame(SOC_dyn).fillna(method = 'bfill')
        SOC_dyn = SOC_dyn.values.flatten()
        SOC_exp_array.append(npy.asarray(SOC_dyn))
        print 'SOC_dyn (with interpolated data)'
        print SOC_dyn

        #change dates to before experiment
        SITE_year0[j] = npy.min(site_T0['Year'])
        date_init = npy.str(npy.int(year_0))
        date_end = npy.str(npy.int(year_end))
        date_init_ss = npy.str(npy.int(year_0 - 30))
        date_end_ss = npy.str(npy.int(year_0 - 1))
	exper_len = npy.int(date_end) - npy.int(date_init) + 1
        SITE_exper_len[j] = exper_len
        clay = npy.mean(site_T0['Clay'])
        SITE_clay[j]=clay
        SITE_date_init[j]=date_init
        SITE_date_end[j]=date_end
        SITE_date_init_ss[j]=date_init_ss
        SITE_date_end_ss[j]=date_end_ss	
        
        #Temperature and moisture variables
	#Steady state
        soil_temp_ss = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_hum_ss  = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
	#Forward
        soil_temp    = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init+"_"+date_end+".txt"
        soil_hum     = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init+"_"+date_end+".txt"
        
        with open(soil_temp_ss) as fileID:
	        # C to K
	        temp_in = npy.array(map(float,fileID))+273.15
	        SITE_temp_in.append(temp_in)

        with open(soil_temp) as fileID:
                # C to K
                temp_t = npy.array(map(float,fileID))+273.15
	        SITE_temp_t.append(temp_t)

        with open(soil_hum_ss) as fileID:
	        #  conversion kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
	        water_in_m3 = npy.array(map(float,fileID))/100
		
        with open(soil_hum) as fileID:
	        #  conversion kg(H2O)/m2(soil) to m3(H2O)/m3(soil)
	        water_t_m3 = npy.array(map(float,fileID))/100
       
        #----------------------------------------------------------------------
        #  determine litter_inc for current site
        #----------------------------------------------------------------------
        # returns ABOVE, BELOW (gC/m2/day) and YEAR of measurements
        SAVE_FILE=-1
        ABOVE,BELOW,YEAR = AB_NanZeroRemover(site_T0,site_T0_name,SAVE_FILE,ROOTDIR)
        #---------------------------------------------------------------------------
        SITE_ABOVE.append(ABOVE)
        SITE_BELOW.append(BELOW)
        ABOVE_mean=npy.mean(ABOVE)
        BELOW_mean=npy.mean(BELOW)
        SITE_ABOVE_mean[j]=ABOVE_mean
        SITE_BELOW_mean[j]=BELOW_mean
        SITE_ERR2_ABOVE_mean[j]=npy.std(ABOVE)**2/len(ABOVE)
        SITE_ERR2_BELOW_mean[j]=npy.std(BELOW)**2/len(BELOW)
        if (SITE_ERR2_ABOVE_mean[j] == 0):
                SITE_ERR2_ABOVE_mean[j]=0.05 # to be checked
                SITE_ERR2_BELOW_mean[j]=0.05
        
        frac_AB_struc=npy.float(frac_array[0]) #fraction of structural on total aboveground litter
        frac_AB_metab=npy.float(frac_array[1]) # fraction of metabolic on total above
        frac_BE_struc=npy.float(frac_array[2]) #fraction of structural on total below
        frac_BE_metab=npy.float(frac_array[3]) #fraction of metabolic on total below


        #mean litter C inputs (gC/m2/day) 
        a_m = ABOVE_mean*frac_AB_metab
        b_m = BELOW_mean*frac_BE_metab
        a_s = ABOVE_mean*frac_AB_struc
        b_s = BELOW_mean*frac_BE_struc

        Err2_am=SITE_ERR2_ABOVE_mean[j]*frac_AB_metab*frac_AB_metab
        Err2_bm=SITE_ERR2_BELOW_mean[j]*frac_BE_metab*frac_BE_metab
        Err2_as=SITE_ERR2_ABOVE_mean[j]*frac_AB_struc*frac_AB_struc
        Err2_bs=SITE_ERR2_BELOW_mean[j]*frac_BE_struc*frac_BE_struc

        litter_inc = npy.array([a_m, b_m, a_s, b_s])    # litter C inputs parameters (gC/m2/day)
        Err2_litter_inc = npy.array([Err2_am,Err2_bm,Err2_as,Err2_bs])
        tot_litter_inc=a_m+b_m+a_s+b_s                  # means total litter carbon inputs per day (gC/m2/day)

        SITE_litterinc[j]      = litter_inc
        SITE_litterinc_err2[j] = Err2_litter_inc

       
        #===================================================
        # SOC,VARIANCE, YEARS with nan and zero removed
        #===================================================
        sc,var,yy0=SOC_NanZeroRemover(site_T0,site_T0_name)
        print sc
        SOC_clean_exp_array.append(sc)
        SOC_clean_exp_variance.append(var)
        SOC_clean_year.append(yy0)
        #===================================================

        # initial pool size; gC/m2
        matrix_cpools = npy.zeros((1,np))
        water_in=water_rh(clay,water_in_m3)
        SITE_water_in.append(water_in)
        SITE_mean_relhum[j]=npy.mean(water_in)
        water_t=water_rh(clay,water_t_m3)
        SITE_water_t.append(water_t)

        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        print 'Used data fort litter_inc from ',len(YEAR),' years over ',len(ABOVE)
        print 'ABOVE_mean: ',ABOVE_mean,'ERR2_above: ',SITE_ERR2_ABOVE_mean[j]
        print 'BELOW_mean: ',BELOW_mean,'ERR2_below: ',SITE_ERR2_BELOW_mean[j]
        print 'BELOW_mean/ABOVE_mean: ',BELOW_mean/ABOVE_mean
        print 'frac_AB_metab,frac_BE_metab,frac_AB_struc,frac_BE_struc: ',frac_AB_metab,frac_BE_metab,
        frac_AB_struc,frac_BE_struc
        print "total litter income (gC/m2/day) ",tot_litter_inc
        print "litter input parameters before 4p1000:", litter_inc 
        print " "
        print "SOC data for fit"
        print " YEAR    SOC     ERROR"
        for k in range(len(sc)):
                print '{0:6d} {1:7.1f} {2:7.1f}'.format(int(yy0[k]),sc[k],math.sqrt(var[k]))
        print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
  
        j+=1

#>>>>>>>> END_OF_INITIALIZATION <<<<<<<<<<<<<<<<<<<<<<<<<<<<        
#exit()

############################################################
#FUNCTIONS
############################################################
#a_matrix
############################################################
def a_matrix(clay):
    a_matrix = npy.zeros((np,np))
    npy.fill_diagonal(a_matrix, -1)
    a_matrix[4,0] = frac_soil_metab_aa                                            # above metabolic to active soil
    a_matrix[4,1] = frac_soil_metab_ab                                            # below metabolic to active soil
    a_matrix[4,2] = frac_soil_struct_aa  * (1- round(lignin_struc_cmatrix[0],2))  # above structural to active soil
    a_matrix[4,3] = frac_soil_struct_ab  * (1 - round(lignin_struc_cmatrix[1],2)) # below structural to active soil
    a_matrix[5,2] = frac_soil_struct_sa * round(lignin_struc_cmatrix[0],2)        # above structural to slow soil
    a_matrix[5,3] = frac_soil_struct_sb * round(lignin_struc_cmatrix[1],2)        # below structural to slow soil
    a_matrix[6,4] = frac_passive_active                                           # active to passive
    a_matrix[5,4] = 1.0 - (0.85-0.68*clay) - a_matrix[6,4]                        # active to slow
    a_matrix[4,5] = frac_active_slow                                              # slow to active
    a_matrix[6,5] = frac_passive_slow                                             # slow to passive
    a_matrix[4,6] = frac_active_passive                                           # passive to active
    a_matrix[5,6] = frac_slow_passive                                             # passive to slow
    a_out=a_matrix
    return a_out

############################################################
#kk_matrix
############################################################
def kk_matrix(tau_array,clay,tsurf_in,tsoil_decomp,litterhum,soilhum_decomp, soilQ10_t):
    kk_matrix = npy.zeros((np,np))
    iabove = 0
    ibelow = 1
    imetabolic = 0
    istructural = 1
    litter_tau=npy.zeros((2))
    litter_tau[imetabolic]  = tau_array[0]
    litter_tau[istructural] = tau_array[1]
    soc_tau = [tau_array[2],tau_array[3],tau_array[4]]
    frozen_respiration_func = 1
    control_temp=npy.zeros((2))
    control_temp[iabove] = control_temp_func(tsurf_in, frozen_respiration_func, soilQ10_t)
    control_temp[ibelow] = control_temp_func(tsoil_decomp, frozen_respiration_func, soilQ10_t)
    control_moist=npy.zeros((2))
    control_moist[iabove] = control_moist_func(litterhum) 
    control_moist[ibelow] = control_moist_func(soilhum_decomp) 
    kk_matrix[0,0] = 1.0/litter_tau[imetabolic]*control_temp[iabove]*control_moist[iabove]
    kk_matrix[1,1] = 1.0/litter_tau[imetabolic]*control_temp[ibelow]*control_moist[ibelow]
    kk_matrix[2,2] = 1.0/litter_tau[istructural]*control_temp[iabove]*control_moist[iabove]*npy.exp(-litter_struct_coef*lignin_struc_cmatrix[0])
    kk_matrix[3,3] = 1.0/litter_tau[istructural]*control_temp[ibelow]*control_moist[ibelow]*npy.exp(-litter_struct_coef*lignin_struc_cmatrix[1])
    kk_matrix[4,4] = 1.0/soc_tau[0]*control_moist[ibelow]*control_temp[ibelow]*(1. - flux_tot_coeff[2]*clay)
    kk_matrix[5,5] = 1.0/soc_tau[1]*control_moist[ibelow]*control_temp[ibelow]
    kk_matrix[6,6] = 1.0/soc_tau[2]*control_moist[ibelow]*control_temp[ibelow]
    return kk_matrix

############################################################
#Spinup
############################################################
def spinup(tau_array,frac_array,clay,temp_in,water_in,soilQ10_t):
    global ABOVE_mean,BELOW_mean, err_above, err_below

    frac_AB_struc=npy.float(frac_array[0]) #fraction of structural on total aboveground litter
    frac_AB_metab=npy.float(frac_array[1]) # fraction of metabolic on total above
    frac_BE_struc=npy.float(frac_array[2]) #fraction of structural on total below
    frac_BE_metab=npy.float(frac_array[3]) #fraction of metabolic on total below


    a_m = ABOVE_mean*frac_AB_metab
    b_m = BELOW_mean*frac_BE_metab
    a_s = ABOVE_mean*frac_AB_struc
    b_s = BELOW_mean*frac_BE_struc

    Err2_am=err_above*frac_AB_metab*frac_AB_metab
    Err2_bm=err_below*frac_BE_metab*frac_BE_metab
    Err2_as=err_above*frac_AB_struc*frac_AB_struc
    Err2_bs=err_below*frac_BE_struc*frac_BE_struc

    litter_inc = npy.array([a_m, b_m, a_s, b_s])    # litter C inputs parameters (gC/m2/day)
    matrix_in_mean=npy.append(litter_inc,[0.,0.,0.])

    for ts in range(0,iforce_recycle):
        tsurf_in = temp_in[ts]
        tsoil_decomp = temp_in[ts]
        litterhum = water_in[ts]
        soilhum_decomp = water_in[ts]
        if (ts == 0):
            kk_ma_mean=kk_matrix(tau_array,clay,tsurf_in,tsoil_decomp,litterhum,soilhum_decomp,soilQ10_t)
        else:
            kk_ma_mean+=kk_matrix(tau_array,clay,tsurf_in,tsoil_decomp,litterhum,soilhum_decomp,soilQ10_t)
    kk_ma_mean=kk_ma_mean/iforce_recycle
    a_ma_mean=a_matrix(clay)
    ss_spinup=-npy.linalg.solve(npy.dot(a_ma_mean,kk_ma_mean),matrix_in_mean)
    return ss_spinup

########################################################################
# forward
########################################################################
def forward(YEARS,init,frac_array,tau_array,clay,temp_t,water_t, soilQ10_t):
    global tsurf_in,tsoil_decomp,litterhum,soilhum_decomp,prior_tau
    global ABOVE_mean, BELOW_mean

    frac_AB_struc=npy.float(frac_array[0]) #fraction of structural on total aboveground litter
    frac_AB_metab=npy.float(frac_array[1]) # fraction of metabolic on total above
    frac_BE_struc=npy.float(frac_array[2]) #fraction of structural on total below
    frac_BE_metab=npy.float(frac_array[3]) #fraction of metabolic on total below

    a_m = ABOVE_mean*frac_AB_metab
    b_m = BELOW_mean*frac_BE_metab
    a_s = ABOVE_mean*frac_AB_struc
    b_s = BELOW_mean*frac_BE_struc

    Err2_am=err_above*frac_AB_metab*frac_AB_metab
    Err2_bm=err_below*frac_BE_metab*frac_BE_metab
    Err2_as=err_above*frac_AB_struc*frac_AB_struc
    Err2_bs=err_below*frac_BE_struc*frac_BE_struc

    litterin = npy.array([a_m, b_m, a_s, b_s])    # litter C inputs parameters (gC/m2/day)

    matrix_next = init                                # starting state vector
    matrix_in=npy.append(litterin,[0.,0.,0.])*dt      # input term I
    a_ma=a_matrix(clay)                               # A matrix
    DY=npy.diff(YEARS)                                # DY (length=#YEARS-1) is the Delta_Time between neighbour years
                                                      # example YEARS=[1990,1993,2000,2018] ==> DY=[3,7,18]
    L=len(DY)
    matrix_cpools_tmean = npy.zeros(L)                # mean total stock of carbon for years form YEARS[1] to YEARS[last]

    #print 'FORWARD'
    #print '-------'
    #print ' '
    #print 'YEARS ',YEARS
    #print 'DY    ',DY
    #print 'INIT  ',matrix_next
    #print 'I     ',matrix_in

    day_start=0
    for x in range(L):               # loop on the number of years after the first.
        DeltaT=DY[x]*one_year
        day_end=day_start+DeltaT
        #print 'Loop on range ',x,' day_start / day_end ',day_start,day_end, ' DY ',DY[x]
        for ts in range(day_start,day_end):  # loop on the # of days corresponding to DY[x] years
            tsurf_in = temp_t[ts]
            tsoil_decomp = temp_t[ts]
    	    litterhum = water_t[ts]
            soilhum_decomp = water_t[ts]
            matrix_current=matrix_next
            kk_ma = kk_matrix(tau_array,clay,tsurf_in,tsoil_decomp,litterhum,soilhum_decomp,soilQ10_t)
            matrix_next = matrix_current + matrix_in + npy.dot(a_ma,npy.dot(kk_ma,matrix_current))*dt
            matrix_cpools_tmean[x]+=npy.sum(matrix_next)
            #print ' ts ',ts, '---------> matrix_cpools_tmean[',x,']=',matrix_cpools_tmean[x]
        day_start=day_end
        matrix_cpools_tmean[x]=matrix_cpools_tmean[x]/DeltaT
        #print 'AVERAGE ',matrix_cpools_tmean[x]

    return matrix_cpools_tmean

############################################################
#control_moist_func
############################################################

def control_moist_func(moist_in):
    moist_coeff=[1.1, 2.4, 0.29]
    moistcont_min=0.25
    moistfunc_result = -moist_coeff[0] * moist_in * moist_in + moist_coeff[1]*moist_in - moist_coeff[2]
    return max(moistcont_min, min(1,moistfunc_result))

############################################################
#control_temp_func
############################################################

#control_temp_plot=open("control_temp.txt","w+")
#temp_plot=open("temp_plot.txt","w+")

def control_temp_func(temp_in, frozen_respiration_func,soilQ10_t):
    soil_Q10 = soilQ10_t[0]	
#   print "SOIL Q10 in control temp",soil_Q10
    tsoil_ref = soilQ10_t[1]
#    print "TEMP ref in control temp",tsoil_ref
    ZeroCelsius = 273.15
    if frozen_respiration_func == 0: #this is the standard ORCHIDEE state
        tempfunc_result= npy.exp(soil_Q10 * (temp_in - (ZeroCelsius+tsoil_ref)) / Q10)
        tempfunc_result= npy.minimum(1.0, tempfunc_result)
    if frozen_respiration_func == 1: #cutoff respiration when T < -1C
        if npy.all(temp_in > ZeroCelsius): #normal as above
            tempfunc_result= npy.exp(soil_Q10 * (temp_in - (ZeroCelsius+tsoil_ref)) / Q10)
        elif npy.all(temp_in > (ZeroCelsius-1.)):
            tempfunc_result = (temp_in-(ZeroCelsius-1.))*npy.exp(soil_Q10*(ZeroCelsius-(ZeroCelsius+tsoil_ref))/Q10)
        else:
            tempfunc_result = 0.0
            tempfunc_result = npy.maximum(npy.minimum(1.0, tempfunc_result), 0)
    if frozen_respiration_func == 2: #cutoff respiration when T < -3C
        if npy.all(temp_in > ZeroCelsius):
            tempfunc_result = npy.exp(soil_Q10 * (temp_in - (ZeroCelsius+tsoil_ref) ) / Q10 )
        elif npy.all(temp_in > (ZeroCelsius - 3.)):
            tempfunc_result = ((temp_in - (ZeroCelsius - 3.))/3.)* npy.exp( soil_Q10 * ( ZeroCelsius - (ZeroCelsius+tsoil_ref) ) / Q10)
        else:
            tempfunc_result = 0.0
    if frozen_respiration_func == 3: #q10 = 100 when below zero
        if npy.all(temp_in > ZeroCelsius):
            tempfunc_result = npy.exp( soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref) ) / Q10)
        else:
            tempfunc_result = npy.exp( 4.605 * ( temp_in - (ZeroCelsius) ) / Q10)* npy.exp( soil_Q10 * ( -tsoil_ref ) / Q10 )
    if frozen_respiration_func == 4: #q10 = 1000 when below zero
        if npy.all(temp_in > ZeroCelsius):
            tempfunc_result = npy.exp(soil_Q10 * ( temp_in - (ZeroCelsius+tsoil_ref) ) / Q10 )
        else:
            tempfunc_result = npy.exp( 6.908 * ( temp_in - (ZeroCelsius) ) / Q10)* npy.exp( soil_Q10 * ( -tsoil_ref ) / Q10)
        
    return  npy.maximum(npy.minimum(1.0, tempfunc_result),0)

        
############################################################
# Q10_opt ######## OBJECTIVE FUNCTION
############################################################
def Jq10_new(q10_new_t):
    global NEW_ITER, clay, temp_in, water_in, prior_tau, SOC_data, SOC_var, YEARS, Current_Site_Index, litter_inc
    global ssp, q10_predict_c_forward, ABOVE_mean, BELOW_mean

    Delta=npy.zeros(len(YEARS))
    j=Current_Site_Index

    ABOVE_mean=SITE_ABOVE_mean[j]
    BELOW_mean=SITE_BELOW_mean[j]

    q10_predict_c = spinup(prior_tau,frac_array,clay,temp_in,water_in,q10_new_t)
    ssp = npy.sum(q10_predict_c)
    q10_predict_c_forward = forward(YEARS,q10_predict_c,frac_array,prior_tau,clay,temp_t,water_t,q10_new_t)	
    Delta[0]=ssp-SOC_data[0]
    Delta[1:]=q10_predict_c_forward-SOC_data[1:]
    Delta2w=Delta**2/SOC_var
    Jq10_new = npy.sum(Delta2w)
    #Jq10_new = npy.sum(Delta**2)
  
    if ( NEW_ITER % CHI2_PRINT_FREQUENCY == 0 ):
        print "==> NEW_ITER ",NEW_ITER," q10_new=",q10_new_t," Jq10_new=",Jq10_new
        print "q10_predict_c_forward",q10_predict_c_forward
        print "SOC_data",SOC_data
        print "Delta",Delta
        print "Delta2w",Delta2w
        print "Error %",npy.sqrt(SOC_var)/SOC_data

    NEW_ITER+=1
    return Jq10_new

#========================================================
#                   MAIN
#========================================================
################################################
# OPTIMIZATION CONSTRAINTS
################################################
#constraints
#def constr_ab(x):
#        return x[0]+x[1]-1
#def constr_be(x):
#        return x[2]+x[3]-1
#con1={'type':'eq','fun':constr_ab}
#con2={'type':'eq','fun':constr_be}
#cons=[con1,con2]


tstart = time.time()

#########
#bounds
#########
bnds=[(0.,2.),(0.,30.)]

Minimizer_trust_constr  = False
Minimizer_SLSQP         = True

q10_sites = npy.zeros((N_sites,2))

# parameter(s) limit(s)
low= None
upp=None


for j in range(N_sites):
        if(Use_Site[j] > 0):
                Current_Site_Index=j
                site       = site_names[j]
                YEARS      = SOC_clean_year[j]
                SOC_data   = SOC_clean_exp_array[j]
                SOC_var    = SOC_clean_exp_variance[j]
                clay       = SITE_clay[j]
                temp_in    = SITE_temp_in[j]
                water_in   = SITE_water_in[j]
                temp_t     = SITE_temp_t[j]
                water_t    = SITE_water_t[j]
		err_above = SITE_ERR2_ABOVE_mean[j]
		err_below = SITE_ERR2_BELOW_mean[j]

                ABOVE_mean = SITE_ABOVE_mean[j]
                BELOW_mean = SITE_BELOW_mean[j]
		
		frac_array=imp_frac[j]
                print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                print '>>>>>>>>>>>>>> Analysis for SITE ',site
                print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                #--------------------------------------------------------------#
                ss=spinup(prior_tau,frac_array,clay,temp_in,water_in,prior_soilQ10_t)
		fwd=forward(YEARS,ss,frac_array,prior_tau,clay,temp_t,water_t,prior_soilQ10_t)
                #--------------------------------------------------------------#
                SITE_SOC_model=npy.sum(ss) #gC/m2
                print 'Stationary solution with priors: '
                print SITE_SOC_model

		SITE_SOC_dyn = npy.concatenate((npy.array([SITE_SOC_model]),fwd))
		print 'SOC dynamics before opt', SITE_SOC_dyn	
		
                print ' '
                if( Minimizer_trust_constr ):
                        #bnds=[(low,upp),(low,upp)]
                        opt_q10_mean=minimize(Jq10_new, frac_array, method='trust-constr', bounds=bnds, options={'disp':True})
                        print "Optimum solution soil_Q10:",opt_q10_mean.x
                        q10_min=opt_q10_mean.x
                        CHI2=Jq10_new(opt_q10_mean.x)
                        print 'TRUST-CONSTR: CHI2 ',CHI2
                elif( Minimizer_SLSQP ):
                        opt_q10_mean=minimize(Jq10_new, prior_soilQ10_t, method='SLSQP', bounds=bnds, options={'disp':True})
                        print "SLSQP: Optimum solution soil_Q10:",opt_q10_mean.x
                        q10_min=opt_q10_mean.x
                        CHI2=Jq10_new(opt_q10_mean.x)
                        print 'SLSQP: CHI2 ',CHI2
			ssp_o = npy.array([ssp])
                        SOC_model = npy.concatenate((ssp_o,q10_predict_c_forward))
			SOC_error=SOC_var
			print 'SOC_error',SOC_error
			XX=npy.stack((YEARS,SOC_model,SOC_data,SOC_error))
			npy.save(out,XX)
        
	                print "SOC model", SOC_model
			print "len SOC model:",len(SOC_model)
			print "len SOC data",len(SOC_data)
                        print "SOC EXPER", SOC_data
                        print "Optimized SOC"
                        print " YEAR    SOCmodel     SOCdata"
                        for k in range(len(SOC_model)-1):
                                print '{0:6d} {1:7.1f} {2:7.1f}'.format(int(YEARS[k]), SOC_model[k],SOC_data[k])
			

                NEW_ITER = 0  
                q10_sites[j] = q10_min

print ' '
for i in range(N_sites):
        if(Use_Site[i] > 0):
                print 'site ',site_names[i],' optimized soil Q10 = ',q10_sites[i]

out2=open("opt_q10Tref_forscript2.10.txt","wb")
npy.save(out2,q10_sites)
out2.close()

optimized_values=q10_sites[~npy.all(q10_sites==0,axis=1)]
opt_out=open("opt_q10_tref2.10.txt","wb")
npy.save(opt_out,optimized_values)
opt_out.close()
out.close()
tend = time.time()
tot_time=tend-tstart
print " "
print " Ok, done. Total time: ",tot_time
