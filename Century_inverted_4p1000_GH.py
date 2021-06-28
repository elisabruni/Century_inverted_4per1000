###########################################################
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

#npy.set_printoptions(threshold=sys.maxsize)
#npy.set_printoptions(linewidth=npy.inf)

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
        BigRelativeError=0.15 # mean percentage variance amongst all sites
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

n_an = 30.

iforce_recycle=30*one_year

#prior_soilQ10 = npy.log(2)
#prior_t = 30.
#prior_soilQ10_t = npy.array([prior_soilQ10,prior_t])

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
ROOTDIR=Root_directory
loc_exp = ROOTDIR+experiment_location

C_input_exp = pd.read_excel(loc_exp)
site_names_all = C_input_exp['ID.Site'].unique()[2:len(C_input_exp)]
site_names_all = map(str, site_names_all)

N_sites_all=len(site_names_all)

#Control plot names
site_T0_array_all=npy.array(['CHNO3_Min', 'COL_T0', 'CREC3_Min', 'FEU_T0', 'JEU2_M0', 'LAJA2_Min', 'LAJA3_Min', 'RHEU1_Min', 'RHEU2_T0','ARAZ_D0_N0', 'ULT_P0_B', 'BROAD_3_Nill', 'FOG_DwN0', 'TREV1_Min','AVRI_T1TR'])



#######
#SITES
#######
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
Use_Site=npy.arange(N_sites_all)

#select sites to be used
Use_Site=npy.zeros(N_sites_all,dtype=npy.int16)
Use_Site[CHNO3]=1
Use_Site[COL]=1
Use_Site[CREC3]=1
Use_Site[FEU]=1
Use_Site[JEU1]=1
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


#Import optimized parameters

##################
#open optimized struc:metab ratios
#################
n=0 # contatore
inp=open('opt_abfractions_forscript2.9.txt','rb')
while inp.read(1):
    inp.seek(-1,1)
    XXmet_all=npy.load(inp)
n+=1

##################
#open optimized Q10 and Tref
#################
n=0 # contatore
inp=open('opt_q10Tref_forscript2.10_corrige.txt','rb')
while inp.read(1):
    inp.seek(-1,1)
    XXQ10_Tref_all=npy.load(inp)
n+=1

##########################
# For non optimized sites (JEU) -> take values from previous site (FEU)
###########################

#non_opt_sites =  npy.where(Use_Site==0)[0] #lines where used_site is 0

non_opt_sites =  npy.where(npy.all(XXmet_all,axis=1)==0)[0] #copy FEU values to JEU
XXmet_all[non_opt_sites]=XXmet_all[non_opt_sites-1]
XXmet=npy.array(XXmet_all)

XXQ10_Tref_all[non_opt_sites]=XXQ10_Tref_all[non_opt_sites-1]
XXQ10_Tref=npy.array(XXQ10_Tref_all)

if(npy.sum(Use_Site)<len(Use_Site)): #if there are unused sites, take them away

	############################
	##Select  met:struc ratios for used only sites
	############################

	XXQ10_Tref_red = npy.array(XXQ10_Tref_all)
	XXmet_red = npy.array(XXmet_all)
	site_names_red = npy.array(site_names_all)
	site_T0_array_red = npy.array(site_T0_array_all)	

	for i,row in enumerate(non_opt_sites):
	        if i==0:
	                del_line = row
			XXQ10_Tref_red = npy.delete(XXQ10_Tref_red,del_line,axis = 0)
	                XXmet_red = npy.delete(XXmet_red,del_line,axis = 0)
			site_names_red = npy.delete(site_names_red,del_line,axis=0)
			site_T0_array_red = npy.delete(site_T0_array_red,del_line,axis=0)
	        else:
	                del_line = row-1
			XXQ10_Tref_red = npy.delete(XXQ10_Tref_red,del_line,axis = 0)
	                XXmet_red = npy.delete(XXmet_red,del_line,axis = 0)
			site_names_red = npy.delete(site_names_red,del_line,axis=0)
			site_T0_array_red = npy.delete(site_T0_array_red,del_line,axis=0)
	
	XXQ10_Tref=npy.array(XXQ10_Tref_red)
	XXmet=npy.array(XXmet_red)
	site_names = npy.array(site_names_red)
	site_T0_array = npy.array(site_T0_array_red)

else:
	XXQ10_Tref= npy.array(XXQ10_Tref_all)
	XXmet = npy.array(XXmet_all)	
	site_names = npy.array(site_names_all)
	site_T0_array = npy.array(site_T0_array_all)


#=========
#Create 4 pools litter income

ns=XXmet.shape[0]
XXstr=1-XXmet

imp_frac = npy.zeros((ns,4))
imp_frac[:,0]=XXstr[:,0]
imp_frac[:,1]=XXmet[:,0]
imp_frac[:,2]=XXstr[:,1]
imp_frac[:,3]=XXmet[:,1]

#========

#==============
print 'USED SITES',site_names
#print imp_frac
#print XXQ10_Tref

#exit()
#=============
#CREATE ARRAYS with data info

#Stationary solution array for each experiment T0
N_sites=len(site_names)

print N_sites

SOC_exp_array=[] #interpolated SOC dynamics experiments
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
SITE_ABOVE_mean=npy.zeros(N_sites)
SITE_BELOW_mean=npy.zeros(N_sites)
SITE_ERR2_ABOVE_mean=npy.zeros(N_sites)
SITE_ERR2_BELOW_mean=npy.zeros(N_sites)
SITE_cov_mean=[]
SITE_mean_relhum=npy.zeros(N_sites)
SITE_litterinc = npy.zeros((N_sites,4))
SITE_litterinc_err2 = npy.zeros((N_sites,4))
SITE_water_in=[]
SITE_temp_in=[]
SITE_water_t=[]
SITE_temp_t=[]
SITE_ABOVE=[]
SITE_BELOW=[]
SITE_TREATMENTS=[]

out_tr=open("SOC_experiments2.txt","wb")

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
	year_30=year_0+30
        missing_years_to30 = npy.int(year_30-year_end-1)

        site_T0_name = site_T0_array[j]
	site_treatments = site_df['ID.Treatment'].unique()[0:len(site_df)]
        site_treatments = map(str, site_treatments)
	#INTERPOLATE each treatment from year0 to year30, fill with Nan
	TREATMENTS=npy.zeros((30,len(site_treatments)))
	count=0
	for i in site_treatments:
                site_T= site_df[(site_df['ID.Treatment'].values == [i])]
                SOC_dyn_T=site_T['SOC']*100 #(gC/m2)
        	if(missing_years_to30>0): #if experiment has less than 30y, fill missing years with Nans
                	empty_ar = npy.empty(missing_years_to30)
                	empty_ar[:]=npy.NaN
			SOC_dyn_T_30=npy.append(SOC_dyn_T,empty_ar)
        	else: #cut experiment to 30th year
                	SOC_dyn_T_30 = npy.array(SOC_dyn_T[0:30])
		TREATMENTS[:,count]=SOC_dyn_T_30
		count+=1

	#TREATMENTS = pd.DataFrame(TREATMENTS,columns = site_treatments)
	SITE_TREATMENTS.append(TREATMENTS)
	npy.save(out_tr,TREATMENTS)
	

        #GET initial years for ss and forward
	site_T0= site_df[(site_df['ID.Treatment'].values == [site_T0_name])]
        SITE_year0[j] = npy.min(site_T0['Year'])
        date_init = npy.str(1980)
        date_end = npy.str(2010)
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
        
        soil_temp_ss = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"temp_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
        soil_hum_ss  = ROOTDIR+'SCRIPT_MODELLI/'+site+'/'+"hum_"+site+"_"+date_init_ss+"_"+date_end_ss+".txt"
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
      
	cov_AB_mean=npy.cov(ABOVE,BELOW)/npy.sqrt(len(ABOVE)) #covariance between ABOVE mean and BELOW mean if no Nans 
	SITE_cov_mean.append(cov_AB_mean)

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

out_tr.close()
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
    frozen_respiration_func = 0
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
def spinup(tau_array,litter_inc,clay,temp_in,water_in,soilQ10_t):
    global ABOVE_mean,BELOW_mean, err_above, err_below

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

############################################################
#Forward
############################################################
def forward(n_an,init,litterin,tau_array,clay,temp,water, soil_Q10):
    global tsurf_t,tsoil_decomp,litterhum,soilhum_decomp,prior_tau # needed!!
    global A_is_constant,dt,one_year # not really needed
    global ABOVE_mean, BELOW_mean

    # flags to determine if A matrix and IN matrix are constant with time or not
    A_is_constant = True
    IN_is_constant = True

    matrix_cpools_tmean = npy.zeros((n_an-1,np))

    length=one_year*n_an
    matrix_cpools_t=npy.zeros((length+1,np))
    matrix_in = npy.zeros(np)
    matrix_cpools_t[0]=init

    for i in range(0,len(litterin)):
        matrix_in[i]=litterin[i]
    if (A_is_constant):
        a_ma=a_matrix(clay)

    for x in range(0,n_an-1):
        matrix_cpools_ymean = npy.zeros(np)
        for ts in range(x*one_year,(one_year*x)+one_year):
            tsurf_t = temp[ts]
            tsoil_decomp = temp[ts]
            litterhum = water[ts]
            soilhum_decomp = water[ts]
            matrix_current= matrix_cpools_t[ts]
            kk_ma = kk_matrix(tau_array,clay,tsurf_t,tsoil_decomp,litterhum,soilhum_decomp, soil_Q10)
            matrix_next = matrix_current + matrix_in + npy.dot(a_ma,npy.dot(kk_ma,matrix_current))*dt
            matrix_cpools_t[ts+1]=matrix_next
            matrix_cpools_ymean += matrix_next
        matrix_cpools_ymean = matrix_cpools_ymean/one_year
	#matrix_cpools_tmean[x] = npy.sum(matrix_cpools_ymean)
        matrix_cpools_tmean[x] = matrix_cpools_ymean
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
# J_new ######## OBJECTIVE FUNCTION
############################################################
def J_new(in_new):

        global NEW_ITER, clay, temp_in, water_in, n_an, Current_Site_Index
        global spinup_c, predict_c, param_opt, ABOVE_mean, BELOW_mean

	predict_c_pools = forward(n_an,spinup_c,in_new,prior_tau,clay,temp_t,water_t,Q10_Tref)
	predict_c = npy.sum(predict_c_pools,axis=1)
	#J_new = abs(target*n_an - sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))
	J_new = abs(target*n_an - npy.sum(predict_c[predict_c.shape[0]-1]-predict_c[0]))

        j=Current_Site_Index
	#ABOVE_mean=SITE_ABOVE_mean[j]
	#BELOW_mean=SITE_BELOW_mean[j]

	NEW_ITER+=1
	if ( NEW_ITER % 100 == 0 ):
	        print "NEW_ITER ",NEW_ITER," in_new=",in_new," J_new=",J_new
	param_opt=in_new
	return J_new


#=============================================
#		MAIN
#=============================================
	
################################################
# OPTIMIZATION CONSTRAINTS
#ordine AM,BM,AS,BS
################################################
def constr1(x):
	global a_constr,ab_ratio
	return x[0]+x[2]-a_constr*(x[1]+x[3])*ab_ratio

def constr2(x):
	global b_constr,ab_ratio
	return b_constr*(x[1]+x[3])*ab_ratio-x[0]-x[2]        


tstart = time.time()

######################################
#Set bounds and constraints for the optimization
#########
bnds=bnds=[(0,10),(0,10),(0,10),(0,10)]

#a_constr=0.8
a_constr=1
#b_constr=1.2
b_constr=1
ab_ratio_var=(1-a_constr)*100
con1={'type':'ineq','fun':constr1}
con2={'type':'ineq','fun':constr2}
cons=[con1,con2]


litterin_sites = npy.zeros((N_sites,4))
SOC_out_all = []
out_mo_pools=open("SOC_model_pools2.txt","wb")
out_mo=open("SOC_model2.txt","wb")
out_lit=open("Litter_income2.txt","wb")
out_priors = open("priors_and_opt_in2.txt","wb")

for j in range(N_sites):
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
	ab_ratio = ABOVE_mean/BELOW_mean
	frac_array=imp_frac[j]
	print 'metabolic:structural fractions',frac_array
	Q10_Tref=XXQ10_Tref[j]
	print 'soilQ10 and Tref', Q10_Tref

	#LITTER INCOME AT SITE	
	#above-below array to calculate uncertainties
	AB_BE_array=npy.array([ABOVE_mean,BELOW_mean])
	
	
	#litter income at site (obs g/m2/day)
        litter_inc = SITE_litterinc[j]
        Err2_litter_inc = SITE_litterinc_err2[j]

	#total litter income at site
	tot_litter_inc = npy.sum(litter_inc)
	tot_err2_litter_inc = npy.sum(Err2_litter_inc)
 
	#to be saved	
	litter_inc_save=npy.append(litter_inc,tot_litter_inc)
	litter_inc_err_save=npy.append(npy.sqrt(Err2_litter_inc),npy.sqrt(tot_err2_litter_inc))

	#litter income prior
	in_opt = npy.array(litter_inc)*(1+0.004)

        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        print '>>>>>>>>>>>>>> Analysis for SITE ',site
        print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
        #--------------------------------------------------------------#
        spinup_c=spinup(prior_tau,litter_inc,clay,temp_in,water_in,Q10_Tref)
	#fwd=forward(n_an,ss,frac_array,prior_tau,clay,temp_t,water_t,Q10_Tref)
        #--------------------------------------------------------------#
        SITE_SOC_model=npy.sum(spinup_c) #gC/m2
        print 'Stationary solution before 4x1000: '
        print SITE_SOC_model
	
	target = SITE_SOC_model*0.004
	n_an = 30
	
	#SITE_SOC_dyn = npy.concatenate((npy.array([SITE_SOC_model]),fwd))
	#print 'SOC dynamics before opt', SITE_SOC_dyn	
	
        print ' '
	print "ABOVE:BELOW ratio allowed to vary by: ", ab_ratio_var, " %"
        opt_mean=minimize(J_new, in_opt, method='SLSQP', bounds=bnds,constraints=cons, options={'disp':True})
	litter_opt = opt_mean.x
        print "SLSQP: Optimum solution:", litter_opt
        total_opt_in=npy.sum(opt_mean.x)

	#optimized litter pools and total litter (save)
	in_opt_save = npy.append(litter_opt,total_opt_in)
	
	#calculate percentage increase/decrease of inputs
        input_change=(total_opt_in-tot_litter_inc)/tot_litter_inc
        print "% change of litter inputs:",input_change

        END=predict_c.shape[0]-1
        C_fin=npy.sum(predict_c[END])
        C_init=npy.sum(predict_c[0])
        SUMO=C_fin-C_init
        Target_reached=(C_fin-C_init)/(C_init*n_an)
        if (Target_reached < 0.005 and Target_reached > 0.003):
                print "Target reached successfully"
                print "Target reached :", Target_reached
        else:
                print "Target not reached"


	############################################################
	#CREATE output file with PREDICTED CARBON over 30 years (standard and 4x100) 
	#       		+ EXPERIMENT SOC filled with Nan for 30 years (T0 + other treatments)
	###########################################################
	predict_c_standard_pools=forward(n_an,spinup_c,litter_inc,prior_tau,clay,temp_t,water_t,Q10_Tref)
	predict_c_standard = npy.sum(predict_c_standard_pools,axis=1)
	predict_c_opt_pools=forward(n_an,spinup_c,litter_opt,prior_tau,clay,temp_t,water_t,Q10_Tref)
	predict_c_opt=npy.sum(predict_c_opt_pools,axis=1)
	
	SOC_model_standard_pools = npy.concatenate((npy.array([spinup_c]),predict_c_standard_pools))
	SOC_model_standard = npy.concatenate((npy.array([SITE_SOC_model]),predict_c_standard))
	SOC_model_opt_pools =npy.concatenate((npy.array([spinup_c]),predict_c_opt_pools))
	SOC_model_opt = npy.concatenate((npy.array([SITE_SOC_model]),predict_c_opt))
	year_out = npy.arange(1,31)
	SOC_pools_out=npy.stack((SOC_model_standard_pools,SOC_model_opt_pools))
	SOC_out=npy.stack((year_out,SOC_model_standard,SOC_model_opt))
	#col_SOC = ['Years','Predicted T0','Predicted 4x1000']
	#SOC_out = pd.DataFrame(npy.transpose(SOC_out),columns = col_SOC)
	SOC_out_all.append(SOC_out)

	npy.save(out_mo_pools,SOC_pools_out)
        npy.save(out_mo,SOC_out)

	#out_c=open("SOC_model_4x1000_"+site_names[j]+".txt","wb")
	#npy.save(out_c,SOC_model_opt)
	#out_c.close()


	############################
	#UNCERTAINTIES
	###########################	
	Uncert_Q = True
        if(Uncert_Q):

                MC_length=50 #set the number of Monte Carlo simulations

                #optimize for n variations of in_opt generatated randomly around the above/below covariance
                opt_parameters_MC=npy.zeros((MC_length,len(litter_inc)))

		#prior above_below 
		ab_be_init_est=AB_BE_array*(1+0.004)

		ABOVE=SITE_ABOVE[j]
		BELOW=SITE_BELOW[j]
		cov_AB_mean=npy.cov(ABOVE,BELOW)/npy.sqrt(len(ABOVE)) #covariance between ABOVE mean and BELOW mean if no Nans
		print 'cov',cov_AB_mean
		if (npy.all(cov_AB_mean)==0): #if covariance is 0, take the mean amongst all sites' covariances to generate cov_AB_mean
			cov_AB_mean=npy.mean(SITE_cov_mean,axis=0)
		print cov_AB_mean

        	frac_AB_struc=npy.float(frac_array[0]) #fraction of structural on total aboveground litter
        	frac_AB_metab=npy.float(frac_array[1]) # fraction of metabolic on total above
        	frac_BE_struc=npy.float(frac_array[2]) #fraction of structural on total below
        	frac_BE_metab=npy.float(frac_array[3]) #fraction of metabolic on total below

		in_rand_param_MC = npy.zeros((MC_length,len(litter_inc)))
		sample_shape=0
                while sample_shape<MC_length: #generate random multinormal until sample_shape=MC_length
                        in_rand_gen=npy.random.multivariate_normal(ab_be_init_est,cov_AB_mean)
			if all(i>0 for i in in_rand_gen): #test if all elements of random array are positive
				#Only positive arrays
				in_rand=in_rand_gen
                        	in_rand_am = in_rand[0]*frac_AB_metab
                        	in_rand_bm = in_rand[1]*frac_BE_metab
                        	in_rand_as = in_rand[0]*frac_AB_struc
                        	in_rand_bs = in_rand[1]*frac_BE_struc
                        	in_rand_param=npy.array([in_rand_am,in_rand_bm,in_rand_as,in_rand_bs])  #array to optimize
				#Save all priors in an array
				in_rand_param_MC[sample_shape]=in_rand_param #add new generated sample to array on rand_in samples
                        	print "inital estimates generated randomly from ab_be_init_est:", in_rand
				#Minimize J_new for the generated array
                        	opt=minimize(J_new,in_rand_param,method='SLSQP',constraints=cons,bounds=bnds,options={'disp':True})
                        	opt_parameters_MC[sample_shape]=opt.x
                        	print "optimum parameter"
                        	print opt_parameters_MC[sample_shape]
				sample_shape+=1
			

                #matrix of the optimum parameters
                print "opt parameters:", opt_parameters_MC

		out_priors_and_opt = npy.stack((in_rand_param_MC,opt_parameters_MC)) #save priors (litterin generated random) and out (optimized)
                npy.save(out_priors,out_priors_and_opt)

                #STANDARD ERROR CALCULATION for the optimized litter inputs to reach 4x1000

                error_AM_opt=npy.std(opt_parameters_MC[:,0])/npy.sqrt(MC_length)
                error_BM_opt=npy.std(opt_parameters_MC[:,1])/npy.sqrt(MC_length)
                error_AS_opt=npy.std(opt_parameters_MC[:,2])/npy.sqrt(MC_length)
                error_BS_opt=npy.std(opt_parameters_MC[:,3])/npy.sqrt(MC_length)
                unc_litter_opt=npy.array([error_AM_opt,error_BM_opt,error_AS_opt,error_BS_opt])
		unc_litter_opt_sum = npy.sum(unc_litter_opt)
                print "Uncertainties:",unc_litter_opt
		
		in_opt_err_save=npy.append(unc_litter_opt,unc_litter_opt_sum)		

		#Error litter = SE per litter in e in opt
		save_lit = npy.stack((litter_inc_save,litter_inc_err_save,in_opt_save,in_opt_err_save))
        	npy.save(out_lit,save_lit)

        NEW_ITER = 0  
        litterin_sites[j] = litter_opt

print ' '
for i in range(N_sites):
        print 'site ',site_names[i],' optimized litter income = ',litterin_sites[i]

optimized_values=litterin_sites[~npy.all(litterin_sites==0,axis=1)]
print 'optimized_values',optimized_values
opt_out=open("opt_litterin2.txt","wb")
npy.save(opt_out,optimized_values)
opt_out.close()
out_mo_pools.close()
out_mo.close()
out_lit.close()
out_priors.close()
tend = time.time()
tot_time=tend-tstart
print " "
print " Ok, done. Total time: ",tot_time
