"""Calculations for the revised charging/collapse logic (manuscript v3, supplement v8, Section S16).
All inputs are stated here; outputs go to model_logic_v8_calcs_record.json."""
import json, math
R=8.314462618; M_S2=0.0641300  # kg/mol
g=3.70; rho=(2000.0,3000.0)
Tb_n, P0 = 717.8, 101325.0      # normal boiling point (CRC; Lide 2007)
Bs=(7750.0,8420.0)               # Clausius-Clapeyron slopes above the boiling point (Section S16)
def Psat(T,B): return P0*math.exp(B*(1/Tb_n-1/T))
def T_at(P,B): return 1/(1/Tb_n - math.log(P/P0)/B)
rec={}
# 1. Gradient of the boiling temperature with depth at the base: dT/dz = T^2/(B z)
h_med=(24.0,28.0); phis=(0.4,0.2,0.1)
grad={}
for phi in phis:
    vals=[]
    for h in h_med:
        d=h/phi
        for r in rho:
            for B in Bs:
                T=T_at(r*g*d,B); vals.append(T*T/(B*d))
    grad[str(phi)]=[round(min(vals),2),round(max(vals),2)]
rec['boiling_T_gradient_K_per_m_at_base']=grad
# 2. Vapor transport vs heat conduction
T=1000.0
vbar=math.sqrt(8*R*T/(math.pi*M_S2))
rec['S2_mean_speed_m_s_1000K']=round(vbar)
L=100.0; kappa=1e-6; yr=3.15576e7
rec['heat_conduction_time_100m_yr']=round(L*L/kappa/yr)
cross={}
for a in (1e-5,1e-3):
    DK=a*vbar/3
    cross[f'aperture_{a*1e6:.0f}um']={'D_K_m2_s':float(f'{DK:.2g}'),
        'crossing_time_100m_days_single_fracture':float(f'{L*L/DK/86400:.2g}'),
        'crossing_time_100m_yr_eps0.01_tau3':float(f'{L*L/(DK*0.01/3)/yr:.2g}')}
rec['vapor_crossing']=cross
# mean free path of S2 at 1 Pa, 1000 K (kinetic diameter ~0.4 nm assumed)
dmol=4e-10; kB=1.380649e-23
rec['S2_mean_free_path_m_1Pa_1000K']=float(f'{kB*T/(math.sqrt(2)*math.pi*dmol**2*1.0):.2g}')
# 3. Sulfur delivery by Knudsen diffusion from FeS (S2 at 0.1-1 Pa) across 100 m
need=100*0.05*2070.0   # kg/m2: 5 vol% native sulfur in 100 m of regolith
rec['sulfur_needed_kg_m2_100m_5volpct']=need
deliv={}
for a in (1e-3,):
    DK=a*vbar/3
    for et in (0.003,0.1):
        for dp in (0.1,1.0):
            J=et*DK*M_S2*dp/(R*T*L)      # kg m-2 s-1
            deliv[f'a{a*1e3:.0f}mm_eps_over_tau{et}_dp{dp}Pa']={'flux_kg_m2_yr':float(f'{J*yr:.2g}'),'time_yr':float(f'{need/(J*yr):.2g}')}
rec['sulfur_delivery_from_FeS']=deliv
# 4. Mobility of native sulfur: vapor pressure at the charging temperatures vs present mean annual temperatures
def P_low(T): return math.exp(24.41-9029.0/T)   # CRC fit, 375-717 K, used below 717.8 K
mob={}
for Tc in (820.0,1050.0):
    for Ts in (270.0,330.0):
        mob[f'{Tc:.0f}K_vs_{Ts:.0f}K']=float(f'{min(Psat(Tc,B) for B in Bs)/P_low(Ts):.2g}')
rec['vapor_pressure_ratio_charging_vs_present']=mob
rec['native_S_vapor_pressure_Pa_500K']=round(P_low(500.0))
# 5. Contact temperature and the latitude test
def Tcontact(Th,Tc,Ih=1.0,Ic=1.0): return (Ih*Th+Ic*Tc)/(Ih+Ic)
rec['contact_T_equal_inertia_K']={f'{Th:.0f}':Tcontact(Th,350.0) for Th in (1300.0,1700.0)}
# if the charging temperature follows half the change in local mean annual surface temperature
Ts=290.0; dTs=21.5
lat={}
for phi,Tc_mid in ((0.4,850.0),(0.2,920.0),(0.1,1000.0)):
    coef=[-B*Ts**2/Tc_mid**2*0.5 for B in Bs]
    ratio=[math.exp(-B*0.5*dTs/Tc_mid**2) for B in Bs]
    lat[str(phi)]={'coef_on_1_over_Ts_K':[round(min(coef)),round(max(coef))],'poleward_depth_ratio':[round(min(ratio),3),round(max(ratio),3)]}
rec['half_weighted_surface_temperature_test']=lat
# 6. Floor-lowering limit by residue: lowering x leaves residue (1-phi)/phi * x; shielding when residue > skin layer
for skin in (0.10,0.22):
    rec[f'max_lowering_before_residue_exceeds_{skin}m']={str(p):round(skin*p/(1-p),3) for p in phis}
# 7. Wang et al. (2020) growth-rate bound per Mercury solar day
rec['wang_vertical_rate_per_solar_day_mm']=round(0.23*175.94/365.25,3)

# 8. Liquid FeS (NIST-JANAF Fe-025, delta-f G in kJ/mol at 1,500 and 1,600 K) and the 5.4 Pa delivery case
dfG={1500:-68.943,1600:-65.382}
rec['FeS_pS2_Pa_JANAF_liquid']={str(T):float(f"{math.exp(2*g*1e3/(R*T))*1e5:.2g}") for T,g in dfG.items()}
DK=1e-3*vbar/3
J=0.1*DK*M_S2*5.4/(R*T*L)
rec['sulfur_delivery_from_FeS']['a1mm_eps_over_tau0.1_dp5.4Pa_1600K']={'flux_kg_m2_yr':float(f'{J*yr:.2g}'),'time_yr':float(f'{need/(J*yr):.2g}')}
# 9. Advance of the sulfur-filled front relative to the heat front (both grow as sqrt(t)):
#    sulfur front L_s^2 ~ 2 D_eff dC t / (rho_S f_S); heat front L_h^2 ~ kappa t (same prefactor convention for both)
rhoS=2070.0; fS=0.05
ratios={}
for et,dp in ((0.1,5.4),(0.1,1.0),(0.003,1.0)):
    Deff=et*DK; dC=dp*M_S2/(R*T)
    ratios[f'eps_over_tau{et}_dp{dp}Pa']=round(math.sqrt(kappa/(2*Deff*dC/(rhoS*fS))),1)
rec['heat_front_over_sulfur_front_distance_ratio_1mm_fS0.05']=ratios
# 10. Heat-crossing times for 60-280 m
rec['heat_crossing_time_yr']={str(x):round(x*x/kappa/yr) for x in (60.0,280.0)}
# 11. Contact temperature if the rock has twice the thermal inertia of the regolith
rec['contact_T_rock_inertia_2x_K']={f'{Th:.0f}':round(Tcontact(Th,350.0,2.0,1.0)) for Th in (1300.0,1700.0)}
# 12. Base-temperature range that spans a twofold range of depth at fixed phi: ln2 = B dT / T^2
rec['dT_for_twofold_depth_K']={f'{Tm:.0f}':[round(math.log(2)*Tm**2/B) for B in Bs] for Tm in (850.0,1000.0)}
# 13. Base temperatures if the confining pressure is the least horizontal stress (0.4-0.6 rho g z)
lhs={}
for phi in phis:
    out=[]
    for h in h_med:
        d=h/phi
        for r in rho:
            for B in Bs:
                Tfull=T_at(r*g*d,B)
                for f in (0.4,0.6):
                    out.append(Tfull-T_at(f*r*g*d,B))
    lhs[str(phi)]=[round(min(out)),round(max(out))]
rec['base_T_lowering_if_confinement_0.4_0.6_rho_g_z_K']=lhs
# 14. Sublimation of native sulfur on sunlit ground (Hertz-Knudsen, S8 vapor, CRC fit below 717.8 K)
M_S8=0.25652
sub={}
for Ts in (300.0,400.0):
    P=P_low(Ts); Jm=P*math.sqrt(M_S8/(2*math.pi*R*Ts)); v=Jm/rhoS
    sub[f'{Ts:.0f}K']={'P_Pa':float(f'{P:.2g}'),'time_to_lose_1mm_s':float(f'{1e-3/v:.2g}')}
rec['sublimation_native_S_vacuum']=sub
# 15. Poleward depth ratio implied by the measured regression coefficient (-291 +/- 99 K) for a 21.5 K drop at 290 K
dinv=1/(290-21.5)-1/290
rec['poleward_ratio_implied_by_coefficient']={'central':round(math.exp(-291*dinv),3),'range_pm_1.96se':[round(math.exp(-(291+1.96*99)*dinv),3),round(math.exp(-(291-1.96*99)*dinv),3)]}

json.dump(rec,open('model_logic_v8_calcs_record.json','w'),indent=1)
print(json.dumps(rec,indent=1))
