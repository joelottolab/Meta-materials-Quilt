# TPV TMM + PV estimate (normal incidence)
import json, argparse, os, numpy as np
h=6.62607015e-34; c=299792458.0; kB=1.380649e-23; qe=1.602176634e-19
def planck_lambda(lam_m,T): a=2*h*c**2/lam_m**5; b=1.0/(np.exp(h*c/(lam_m*kB*T))-1.0); return a*b
def nk_const(n,k=0.0): return lambda lam_um:(n,k)
DEFAULT_MATERIALS={"Air":nk_const(1,0),"SiO2":nk_const(1.45,0),"TiO2":nk_const(2.4,0.02),
                   "SiC":nk_const(2.6,0.05),"W":nk_const(3.5,3.0),"BackReflector":nk_const(0.2,7.0)}
def load_nk(nm): 
    if nm in DEFAULT_MATERIALS: return DEFAULT_MATERIALS[nm]
    raise ValueError("Unknown "+nm)
def tmm_R_T_normal(lam, layers):
    lam=np.array(lam); k0=2*np.pi/lam
    n0=layers[0][0]; ns=layers[-1][0]
    M11=np.ones_like(lam,dtype=complex); M12=M21=np.zeros_like(lam,dtype=complex); M22=np.ones_like(lam,dtype=complex)
    for j in range(1,len(layers)-1):
        nj=layers[j][0]; dj=layers[j][1]; beta=k0*nj*dj; cj=np.cos(beta); sj=1j*np.sin(beta); Yj=nj
        M11,M12,M21,M22=(M11*cj+(sj/Yj)*M12, M11*sj*Yj+M12*cj, M21*cj+(sj/Yj)*M22, M21*sj*Yj+M22*cj)
    eta0=n0; etas=ns; Yin=(M11*etas+M12)/(M21*etas+M22); r=(Yin-eta0)/(Yin+eta0); t=2*eta0/(Yin+eta0)
    R=np.real(r*np.conj(r)); T=np.real((etas/eta0)*t*np.conj(t)); R=np.clip(R,0,1); T=np.clip(T,0,1); A=np.clip(1-R-T,0,1); return R,T,A
def build_layers(desc,lam):
    layers=[]; 
    for L in desc:
        nk=load_nk(L["material"]); n_arr=np.array([nk(x)[0]+1j*nk(x)[1] for x in lam]); layers.append([n_arr,float(L["thickness_um"])])
    return layers
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); args=ap.parse_args()
    cfg=json.load(open(args.config)); out=os.path.dirname(args.config); os.makedirs(out, exist_ok=True)
    lam=np.linspace(cfg["simulation"]["lambda_min_um"],cfg["simulation"]["lambda_max_um"],cfg["simulation"]["points"]); T=cfg["simulation"]["emitter_temperature_K"]
    Re,Te,Ae=tmm_R_T_normal(lam, build_layers(cfg["emitter_stack"], lam)); emiss=Ae
    Rf,Tf,Af=tmm_R_T_normal(lam, build_layers(cfg["filter_stack"], lam))
    lam_m=lam*1e-6; B=planck_lambda(lam_m,T); spectral=np.pi*emiss*B*1e-6; post=spectral*Tf
    Eg=cfg["pv_cell"]["bandgap_eV"]; lam_g=1.239841984/Eg; E= (h*c)/(lam_m); photons=post/(E*1e6)
    EQE=np.where(lam<=lam_g, cfg["pv_cell"]["external_quantum_efficiency"], 0.0)
    Jsc=qe*np.trapz(photons*EQE, lam); Voc=cfg["pv_cell"]["voltage_fraction_Voc"]*Eg
    P=Jsc*Voc; inc=np.trapz(spectral, lam); inc_post=np.trapz(post, lam)
    json.dump({"Jsc_A_per_m2":float(Jsc),"Voc_V":float(Voc),"PowerDensity_W_per_m2":float(P),
               "OverallEff_vs_Emitter":float(P/inc) if inc>0 else 0.0,
               "FilterPassEff":float(P/inc_post) if inc_post>0 else 0.0,
               "IncPower_Emitter_W_per_m2":float(inc),"IncPower_PostFilter_W_per_m2":float(inc_post),
               "Lambda_g_um":float(lam_g)}, open(os.path.join(out,"tpv_summary.json"),"w"), indent=2)
    np.savetxt(os.path.join(out,"tpv_spectra.csv"), np.column_stack([lam, emiss, Tf, post]), delimiter=",",
               header="lambda_um,EmitterEmissivity,FilterT,SpectralPowerPostFilter_W_per_m2_um", comments="")
    print("TPV done:", os.path.join(out,"tpv_spectra.csv"))
if __name__=="__main__": main()
