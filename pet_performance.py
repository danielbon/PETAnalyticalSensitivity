#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# evaluate fitting with just one point (peak sensitivity)
# include off axis sensitivity
from __future__ import print_function
from __future__ import annotations
import numpy as np
from numpy import log as ln
from statistics import mean
import math
from array import array
from lmfit import Model
from nistcalculators.xcom import xcom
from scipy.constants import e, h, hbar, alpha, c, m_e

# Several constants factored into a single variable.
f = (hbar * alpha / m_e / c)**2 / 2

class photopeak_efficiency:

    def __init__(self, scintillator, scint_length, energy_inc = 511e3): #scint_length in mm and energy_inc in eV
        self.scintillator = scintillator
        self.scint_length = scint_length
        self.energy_inc = energy_inc
        self.calculate_probabilities()

    def calculate_probabilities(self):
        m_total = 0.
        area_total = 0.
        for i in range(len(self.scintillator)):
            m_total += self.scintillator[i][1]*self.scintillator[i][2]*self.scint_length
            area_total += self.scintillator[i][2]
        formulas_scint=[]
        weights_scint=[]
        for i in range(len(self.scintillator)): 
            weights_scint.append((self.scintillator[i][1]*self.scintillator[i][2]*self.scint_length)/m_total)
            formulas_scint.append(self.scintillator[i][0])
        self.scint_material = xcom.MaterialFactory.mix_formulas(formulas=formulas_scint, weights=weights_scint)
        self.eff_density = m_total/(area_total*self.scint_length)
        self.mass_att_coefficients = xcom.calculate_attenuation(self.scint_material, [self.energy_inc])
        mu_cohe_inc = self.mass_att_coefficients["coherent"][0]*self.eff_density
        mu_comp_inc = self.mass_att_coefficients["incoherent"][0]*self.eff_density
        mu_phe_inc = self.mass_att_coefficients["photoelectric"][0]*self.eff_density
        mu_total_inc = self.mass_att_coefficients["total"][0]*self.eff_density # coherent scattering does not deposit energy
        #mu_total_inc = self.mass_att_coefficients["total_without_coherent"][0]*self.eff_density # coherent scattering does not deposit energy
        #self.interact_diam = 2.*math.sqrt(area_total/np.pi)*0.1 #given an area in mm, convert into an interaction radius in cm
        self.interact_diam = math.sqrt(area_total)*0.1 #given an area in mm, convert into an interaction radius in cm
        print("self.interact_diam: "+ str(self.interact_diam))
        self.length = self.scint_length*0.1 #change from mm to cm
        dL = 0.1 # step of 0.1 cm
        n_layers = int(self.length/dL)
        n_deg = 90 # step of 2 degrees over 180 degrees (1 pi)
        d_deg = np.pi/n_deg
        self.P_no_int = math.exp(-mu_total_inc*self.length) # probability of no interaction
        self.P_phe = (1 - math.exp(-mu_total_inc*self.length))*mu_phe_inc/mu_total_inc # probability of first interaction is photoelectric absorption
        P_comp_esc = 0. # number of photons that escape the detector after a compton scattering
        P_comp_phe = 0. # number of photoelectric interactions after a compton scattering
        P_comp_comp = 0. # number of secondary compton interactions after a first compton scattering
        P_comp = 0.
        P_comp_sum = 0.
        P_comp_check = (1. - math.exp(-mu_total_inc*self.length))*mu_comp_inc/mu_total_inc
        print("P_comp_check: "+ str(P_comp_check))
        for i in range(0, n_layers, 1):
            DOI = i*dL
            print("DOI: "+ str(DOI))
            N_int_DOIdL = 1. - math.exp(-mu_total_inc*(DOI+dL)) # number of interactions from 0 to DOI+dL
            N_int_DOI = 1. - math.exp(-mu_total_inc*(DOI))  # number of interactions in from 0 to DOI
            P_comp = (N_int_DOIdL-N_int_DOI)*mu_comp_inc/mu_total_inc # number of Compton interactions from DOI to DOI+dL
            P_comp_sum += P_comp
            P_esc2 = 0. 
            P_phe2 = 0. 
            P_comp2 = 0. 
            Sum_KN_dsigma_dTheta = 0.
            for j in range(0, n_deg-1, 1):
                degree = float(j*d_deg)
                #print("degree: "+ str(degree))
                W = self.get_KN_W(degree, self.energy_inc)
                E_sct = self.energy_inc*W
                #print("E_sct: "+ str(E_sct))
                length_sct = self.get_length_sct(degree, DOI+dL/2., 1e-8)
                #print("length_sct: "+ str(length_sct))
                self.mass_att_coefficients = xcom.calculate_attenuation(self.scint_material, [E_sct])
                mu_cohe_sct = self.mass_att_coefficients["coherent"][0]*self.eff_density
                mu_comp_sct = self.mass_att_coefficients["incoherent"][0]*self.eff_density
                mu_phe_sct = self.mass_att_coefficients["photoelectric"][0]*self.eff_density
                mu_total_sct = self.mass_att_coefficients["total"][0]*self.eff_density
                #mu_total_sct = self.mass_att_coefficients["total_without_coherent"][0]*self.eff_density
                Sum_KN_dsigma_dTheta += self.get_KN_dsigma_dTheta(degree, self.energy_inc)
                P_esc2 += self.get_KN_dsigma_dTheta(degree, self.energy_inc)*math.exp(-mu_total_sct*length_sct)# photons that escape the detector
                P_phe2 += self.get_KN_dsigma_dTheta(degree, self.energy_inc)*(1-math.exp(-mu_total_sct*length_sct))*mu_phe_sct/mu_total_sct # number of photoelectric innteractions
                P_comp2 += self.get_KN_dsigma_dTheta(degree, self.energy_inc)*(1-math.exp(-mu_total_sct*length_sct))*mu_comp_sct/mu_total_sct # number of compton interactions
            P_esc2 = P_esc2/Sum_KN_dsigma_dTheta
            P_phe2 = P_phe2/Sum_KN_dsigma_dTheta
            P_comp2 = P_comp2/Sum_KN_dsigma_dTheta
            P_comp_esc += P_comp*P_esc2
            P_comp_phe += P_comp*P_phe2
            P_comp_comp += P_comp*P_comp2
        print("P_comp_sum: "+ str(P_comp_sum))
        print("P_comp_check: "+ str(P_comp_check))
        print("P_phe: "+ str(self.P_phe))
        print("P_n_int: "+ str(1.-self.P_no_int)) #number of interactions inside the detector
        print("P_comp_esc: "+ str(P_comp_esc))
        print("P_comp_phe: "+ str(P_comp_phe))
        print("P_comp_comp: "+ str(P_comp_comp))
        print("P_no_int: "+ str(self.P_no_int))
        self.P_upper_bound = 1 - self.P_no_int - P_comp_esc
        self.P_lower_bound = self.P_phe + P_comp_phe
        self.P_mid_bound = self.P_phe + P_comp_phe + P_comp_comp


    def get_length_sct(self, theta, DOI, tolerance=1e-8): #DOI is the depth of interation
        l_sct = 0.  # Default value
        a = self.interact_diam #0.3
        theta1 = np.arctan(a/2./(self.length - DOI))
        theta2 = np.pi - np.arctan(a/2./DOI)
        if 0. <= theta <= theta1:
            l_sct = (self.length - DOI)/np.cos(theta)
        elif theta1 < theta <= np.pi/2.:
            l_sct = (a/2.)/np.sin(theta)
        elif np.pi/2 < theta < theta2:
            l_sct = (a/2.)/np.sin(np.pi-theta)
        elif theta2 <= theta <= np.pi:
            l_sct = DOI/np.cos(np.pi-theta)
        return l_sct

    def get_length_sct_infinite_sides(self, theta, DOI, tolerance=1e-8): #DOI is the depth of interation
        l_sct = 0.  # Default value
        if 0. <= theta < np.pi/2:
            l_sct = (self.length - DOI) / np.cos(theta)
        elif np.pi/2 < theta <= np.pi:
            l_sct = DOI / np.cos(np.pi-theta)
        elif abs(theta - np.pi/2) < tolerance:
            l_sct = math.inf
        return l_sct

    def get_KN_W(self, theta, energy_i = [511e3]): # KN = Klein Nishina; energy is given in eV  
        # Incoming photon frequency (s-1) and wavelength (m).
        nu = energy_i * e / h
        wav_i = c / nu
        # Scattered photon wavelength (m).
        wav_s = wav_i + h / m_e / c * (1 - np.cos(theta))
        W = wav_i / wav_s
        return W # equal to E_incoming/E_scattered

    def get_KN_dsigma_dOmega(self, theta, energy = [511e3]):
        W = self.get_KN_W(theta,energy)
        # Differential cross section given by the Klein-Nishina formula.
        dsigma_dOmega = f * W**2 * (W + 1/W - np.sin(theta)**2)
        return dsigma_dOmega

    def get_KN_dsigma_dTheta(self, theta, energy = [511e3]):
        W = self.get_KN_W(theta,energy)
        # Differential cross section given by the Klein-Nishina formula.
        dsigma_dTheta = f * W**2 * (W + 1/W - np.sin(theta)**2)*2*np.pi*np.sin(theta)
        #dsigma_dTheta=(1/(2-cos(x)))^2*(1/(2-cos(x))+2-cos(x)-sin(x)^2)*sin(x)
        return dsigma_dTheta

    def get_attenuation_coefficient(self, interaction="total", energy = [511e3]):
        attenuation_coefficient = xcom.calculate_attenuation(self.scint_material, energy)
        return attenuation_coefficient[interaction]*self.eff_density

    # pde stands for photopeak detection efficiency
    def get_pde(self, phi=1, interaction="total", energy=[511e3]): #phi=1 when no energy window selection is applied
        pde = 0.
        if interaction == "lower_bound":
            pde = self.P_lower_bound
        elif interaction == "upper_bound":
            pde = self.P_upper_bound        
        elif interaction == "photoelectric":
            pde = self.P_phe
        elif interaction == "mid_bound":
            pde = self.P_mid_bound
        elif interaction == "total": 
            pde = phi*(1 - math.exp(-self.get_attenuation_coefficient(interaction, energy)*self.length))
        else: 
            print("Please choose correct answer. Total being used.")
            pde = phi*(1 - math.exp(-self.get_attenuation_coefficient(interaction, energy)*self.length))
        return pde

    # get the efficiency of an energy window with limits [Emin, Emax] for a gaussian photopeak distribution with mean energy Emen and energy resolution Eres
    def get_ew_efficiency(Eres, Emean, Emin, Emax):
        fwhm = Emean*Eres
        sigma = fwhm / (2.355)  # Convert FWHM of the energy resolution to sigma
        p = 0.5 * (math.erf((Emax - Emean) / sigma / (2 ** 0.5)) - math.erf((Emin - Emean) / sigma / (2 ** 0.5)))
        return p
        

class solid_angle:
    def __init__(self, ring_diam, axial_ext, n_rings, n_mod_per_ring, n_el_per_mod, entrance_area_el):
        self.n_rings = n_rings
        self.n_mod_per_ring = n_mod_per_ring
        self.n_el_per_mod = n_el_per_mod
        self.ring_diam = ring_diam
        self.axial_ext = axial_ext
        self.entrance_area_el = entrance_area_el
        self.cylinder_area = 2*math.pi*(self.ring_diam/2)*self.axial_ext

    def get_solid_angle_center(self):
        return self.axial_ext/math.sqrt(self.axial_ext**2+4*(self.ring_diam/2)**2)

    # from L. Eriksson et al., “An investigation of sensitivity limits in PET scanners,” NIM A, vol. 580, no. 2, pp. 836–842, 2007.
    def get_solid_angle_correlated_peak_ref(self):
        OmegaAxCorRef = self.axial_ext/np.sqrt((self.axial_ext/2.)**2+(self.ring_diam/2)**2) 
        return OmegaAxCorRef

    def get_solid_angle_uncorrelated(self, axial_pos):
        OmegaAxUncor = 0.5*(((self.axial_ext/2.)-axial_pos)/math.sqrt(((self.axial_ext/2.)-axial_pos)**2+(self.ring_diam/2.)**2)+((self.axial_ext/2.)+axial_pos)/math.sqrt(((self.axial_ext/2.)+axial_pos)**2+(self.ring_diam/2.)**2))          
        return OmegaAxUncor

    def get_solid_angle_correlated(self, axial_pos):
        OmegaAxCor = ((self.axial_ext/2.)-np.abs(axial_pos))/np.sqrt(((self.axial_ext/2.)-np.abs(axial_pos))**2+(self.ring_diam/2.)**2)        
        return OmegaAxCor

    def get_solid_angle3D_uncorrelated(self, axial_pos, radial_pos):
        OmegaAxUncor = 0.5*(((self.axial_ext/2.)-axial_pos)/math.sqrt(((self.axial_ext/2.)-axial_pos)**2+(self.ring_diam/2.-radial_pos)**2)+((self.axial_ext/2.)+axial_pos)/math.sqrt(((self.axial_ext/2.)+axial_pos)**2+(self.ring_diam/2.-radial_pos)**2))          
        return OmegaAxUncor

    def get_solid_angle3D_correlated(self, axial_pos, radial_pos):
        OmegaAxCor = ((self.axial_ext/2.)-np.abs(axial_pos))/np.sqrt(((self.axial_ext/2.)-np.abs(axial_pos))**2+(self.ring_diam/2.-radial_pos)**2)        
        return OmegaAxCor


class sensitivity:

    def __init__(self, scanner_name, scintillator, scint_length, ring_diam, axial_ext, n_rings, n_mod_per_ring, n_el_per_mod, entrance_area_el):
        self.scanner_name = scanner_name
        self.scintillator = scintillator
        self.scint_length = scint_length
        self.n_rings = n_rings
        self.n_mod_per_ring = n_mod_per_ring
        self.n_el_per_mod = n_el_per_mod
        self.ring_diam = ring_diam
        self.axial_ext = axial_ext
        self.entrance_area_el = entrance_area_el
        self.cylinder_area = 2*math.pi*(self.ring_diam/2)*self.axial_ext
        self.packing_frac = self._calculate_packing_frac()
        self.packing_frac_ref = self._calculate_packing_frac_ref()
        self.scint_pde = photopeak_efficiency(scintillator, scint_length)
        self.scint_sa = solid_angle(ring_diam, axial_ext, n_rings, n_mod_per_ring, n_el_per_mod, entrance_area_el)

    def __str__(self):
        return f"{self.scanner_name} uses {self.scintillator} scintillator with {self.scint_length} mm"

    def get_cylinder_area(self):
        return self.cylinder_area

    def _calculate_packing_frac(self):
        return self.n_rings * self.n_mod_per_ring * self.n_el_per_mod * self.entrance_area_el / self.cylinder_area

    def _calculate_packing_frac_ref(self):
        return self.n_rings*self.n_mod_per_ring*self.n_el_per_mod*self.entrance_area_el*self.scint_length/(self.axial_ext*math.pi*((self.ring_diam/2+self.scint_length)**2-(self.ring_diam/2)**2))

    def get_packing_frac(self):
        return self.packing_frac

    def force_packing_frac(self, PF):
        self.packing_frac = PF 

    # from L. Eriksson et al., “An investigation of sensitivity limits in PET scanners,” NIM A, vol. 580, no. 2, pp. 836–842, 2007.
    def get_packing_frac_ref(self):
        return self.packing_frac_ref

    # from L. Eriksson et al., “An investigation of sensitivity limits in PET scanners,” NIM A, vol. 580, no. 2, pp. 836–842, 2007.
    def get_peak_sensitivity_ref(self): #pure beta plus emission (two gammas back-to-back)
        return 100*self.scint_sa.get_solid_angle_correlated_peak_ref()*self.get_packing_frac_ref()**2*self.scint_pde.get_pde(1, "total", [511e3])**2

    def get_sensitivity3D_correlated(self, axial_pos, radial_pos, phi, interaction, energy1, energy2): #pure beta plus emission (two gammas back-to-back)
        return 100*self.scint_sa.get_solid_angle3D_correlated(axial_pos, radial_pos)*self.get_packing_frac()*self.scint_pde.get_pde(phi, interaction,energy1)*self.scint_pde.get_pde(phi, interaction,energy2)

    def get_sensitivity3D_uncorrelated(self, axial_pos, radial_pos, phi, interaction, energy1, energy2): #e.g. one annihilation 511 keV photon and a gamma emmission within time window coincidence
        return 100*(self.scint_sa.get_solid_angle3D_uncorrelated(axial_pos, radial_pos)*self.get_packing_frac())**2*self.scint_pde.get_pde(phi, interaction, energy1)*self.scint_pde.get_pde(phi, interaction, energy2)
        #return 100*(self.scint_sa.get_solid_angle_uncorrelated(axial_pos, radial_pos)*self.get_packing_frac())*self.scint_pde.get_pde(phi, interaction, energy1)*self.scint_pde.get_pde(phi, interaction, energy2)

    def compute_axial3D(self, interaction, emitter, axial_pos, radial_pos, phi=1):  #phi=1 when no energy window selection is applied
        if emitter == "sodium-22": 
            return self.sensitivity3D_sodium_22(axial_pos, radial_pos, phi, interaction)
        elif emitter == "511 keV": 
            return self.sensitivity3D_511keV(axial_pos, radial_pos, phi, interaction)
        else: 
            return self.sensitivity3D_511keV(axial_pos, radial_pos, phi, interaction)

    def compute_axial3D_array(self, interaction, emitter, axial_posis, radial_pos, phi=1):  #phi=1 when no energy window selection is applied
        output = []
        for axial_posi in axial_posis:
            output.append(self.compute_axial3D(interaction, emitter, axial_posi, radial_pos, phi))
        return output

    def sensitivity3D_sodium_22(self, axial_posis, radial_pos, phi, interaction = "total"):
        energy_annih = [511e3]
        energy_gamma = [1274e3]
        branch_frac = 1. #0.9 is the branching fraction of sodium-22
        return branch_frac*(self.get_sensitivity3D_correlated(axial_posis, radial_pos, phi, interaction, energy_annih, energy_annih)+2*self.get_sensitivity3D_uncorrelated(axial_posis, radial_pos, phi, interaction, energy_annih, energy_gamma))

    def sensitivity3D_511keV(self, axial_posis, radial_pos, phi, interaction = "total"):
        energy_annih = [511e3]
        branch_frac = 1.0
        return branch_frac*(self.get_sensitivity3D_correlated(axial_posis, radial_pos, phi, interaction, energy_annih, energy_annih))

    def get_sensitivity_correlated(self, axial_pos, phi, interaction, energy1, energy2): #pure beta plus emission (two gammas back-to-back)
        return 100*self.scint_sa.get_solid_angle_correlated(axial_pos)*self.get_packing_frac()*self.scint_pde.get_pde(phi, interaction,energy1)*self.scint_pde.get_pde(phi, interaction,energy2)

    def get_sensitivity_uncorrelated(self, axial_pos, phi, interaction, energy1, energy2): #e.g. one annihilation 511 keV photon and a gamma emmission within time window coincidence
        return 100*(self.scint_sa.get_solid_angle_uncorrelated(axial_pos)*self.get_packing_frac())**2*self.scint_pde.get_pde(phi, interaction, energy1)*self.scint_pde.get_pde(phi, interaction, energy2)
        #return 100*(self.scint_sa.get_solid_angle_uncorrelated(axial_pos)*self.get_packing_frac())*self.scint_pde.get_pde(phi, interaction, energy1)*self.scint_pde.get_pde(phi, interaction, energy2)
    
    def compute_axial(self, interaction, emitter, axial_pos, phi=1):  #phi=1 when no energy window selection is applied
        if emitter == "sodium-22": 
            return self.sensitivity_sodium_22(axial_pos, phi, interaction)
        elif emitter == "511 keV": 
            return self.sensitivity_511keV(axial_pos, phi, interaction)
        else: 
            return self.sensitivity_511keV(axial_pos, phi, interaction)

    def compute_axial_array(self, interaction, emitter, axial_posis, phi=1):  #phi=1 when no energy window selection is applied
        output = []
        for axial_posi in axial_posis:
            output.append(self.compute_axial(interaction, emitter, axial_posi, phi))
        return output

    def sensitivity_sodium_22(self, axial_posis, phi, interaction = "total"):
        energy_annih = [511e3]
        energy_gamma = [1274e3]
        branch_frac = 1. #0.9 is the branching fraction of sodium-22
        return branch_frac*(self.get_sensitivity_correlated(axial_posis, phi, interaction, energy_annih, energy_annih)+2*self.get_sensitivity_uncorrelated(axial_posis, phi, interaction, energy_annih, energy_gamma))

    def sensitivity_511keV(self, axial_posis, phi, interaction = "total"):
        energy_annih = [511e3]
        branch_frac = 1.0
        return branch_frac*(self.get_sensitivity_correlated(axial_posis, phi, interaction, energy_annih, energy_annih))

    def fit_axial_curve(self, interaction, emitter, axial_posis, ref_data, ref_data_err):
        if emitter == "sodium-22":
            gmodel = Model(self.sensitivity_sodium_22)#fit with default interaction = "total"
            fit_params = gmodel.make_params(phi=0.01)
            fit_params['phi'].set(min=0.001, max=0.999)
            result = gmodel.fit(np.array(ref_data), fit_params, axial_posis=axial_posis, weights=1.0/np.array(ref_data_err))
            #print(result.fit_report())
            phi = result.params['phi'].value
            phi_err = result.params['phi'].stderr
            return phi, phi_err, result.redchi
        else:
            gmodel = Model(self.sensitivity_511keV) #fit with default interaction = "total"
            fit_params = gmodel.make_params(phi=0.01)
            fit_params['phi'].set(min=0.001, max=0.999)
            result = gmodel.fit(np.array(ref_data), fit_params, axial_posis=axial_posis, weights=1.0/np.array(ref_data_err))
            #print(result.fit_report())
            phi = result.params['phi'].value
            phi_err = result.params['phi'].stderr
            return phi, phi_err, result.redchi


