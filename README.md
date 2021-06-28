# Century_inverted_4per1000
Daily Century SOC model with inversion algorithm to calculate carbon inputs required to increase SOC stocks by 4 per 1000.
Code created by Elisa Bruni:
The original Century model was developed by Parton et al. (1988).
The matricial version of the model was provided by Yuanyuan Huang (Yuanyuan Huang yuanyuanhuang2011@gmail.com) in matlab code.
Translation to python and further code developments (i.e. parameters calibration, inversion algorithm and uncertainties quantification) were made by Elisa Bruni (ebruni93@gmail.com).

Code used to generate simulations in : Bruni et al. (2021) https://doi.org/10.5194/bg-18-1-2021

Simulations require site level long-term observed data.

The repository includes:
- two calibration files: 
    Century_calibration_ab1_GH.py, for metabolic : structural ratios calibration  
    Century_calibration_qT2_GH.py, for Q10 and reference temperature calibration
- one main file
    Century_inverted_4p1000_GH.py, with the inversion algorithm to calculate carbon inputs required to increase SOC stocks by 4 per 1000 per year. 
    It includes uncertainties quantification.

For questions, comments, or inquiries, please contact Elisa Bruni ebruni93@gmail.com
