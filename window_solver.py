import json, argparse, os, numpy as np
def nk_const(n,k=0.0): return lambda lam:(n,k)
def build_materials(base,state):
    mats={k:nk_const(v["n"],v["k"]) for k,v in base.items()}
    mats["WO3"]=nk_const(state["WO3"]["n"],state["WO3"]["k"])
    mats["NiO"]=nk_const(state["NiO"]["n"],state["NiO"]["k"])
    return mats
def tmm_RT(lam,layers,mats):
    lam=np.array(lam); k0=2*np.pi/lam
    def n_of(name): return np.array([mats[name](L)[0]+1j*mats[name](L)[1] for L in lam])
    n0=n_of(layers[0]["material"]); ns=n_of(layers[-1]["material"])
    M11=np.ones_like(lam,dtype=complex); M12=M21=np.zeros_like(lam,dtype=complex); M22=np.ones_like(lam,dtype=complex)
    for j in range(1,len(layers)-1):
        nj=n_of(layers[j]["material"]); dj=layers[j]["thickness_um"]; beta=k0*nj*dj; cj=np.cos(beta); sj=1j*np.sin(beta); Yj=nj
        M11,M12,M21,M22=(M11*cj+(sj/Yj)*M12, M11*sj*Yj+M12*cj, M21*cj+(sj/Yj)*M22, M21*sj*Yj+M22*cj)
    eta0=n0; etas=ns; Yin=(M11*etas+M12)/(M21*etas+M22); r=(Yin-eta0)/(Yin+eta0); t=2*eta0/(Yin+eta0)
    R=np.real(r*np.conj(r)); T=np.real((etas/eta0)*t*np.conj(t)); R=np.clip(R,0,1); T=np.clip(T,0,1); A=np.clip(1-R-T,0,1); return R,T,A
def luminous_transmittance(lam,T): V=np.exp(-0.5*((lam-0.555)/0.1)**2); return np.trapz(T*V,lam)/np.trapz(V,lam)
def simulate_kinetics(profile,tau_on,tau_off,dt=0.05):
    times=[profile[0]["t"]]; Vs=[profile[0]["V"]]
    for i in range(1,len(profile)):
        t0=profile[i-1]["t"]; v0=profile[i-1]["V"]; t1=profile[i]["t"]; v1=profile[i]["V"]
        steps=max(1,int((t1-t0)/dt))
        for s in range(1,steps+1):
            times.append(t0+s*dt); Vs.append(v0+(v1-v0)*s/steps)
    times=np.array(times); Vs=np.array(Vs); X=np.zeros_like(times); X[0]=0.0
    for i in range(1,len(times)):
        dt_i=times[i]-times[i-1]; target=1.0 if Vs[i]>0.5 else (0.0 if Vs[i]<-0.5 else X[i-1])
        tau=tau_on if target>X[i-1] else tau_off; X[i]=X[i-1]+(target-X[i-1])*(1-np.exp(-dt_i/max(tau,1e-3)))
    return times,Vs,X
def main():
    ap=argparse.ArgumentParser(); ap.add_argument("--config",required=True); args=ap.parse_args()
    cfg=json.load(open(args.config)); out=os.path.dirname(args.config); os.makedirs(out, exist_ok=True)
    lam=np.linspace(cfg["simulation"]["lambda_min_um"],cfg["simulation"]["lambda_max_um"],cfg["simulation"]["points"])
    layers=cfg["stack"]; base=cfg["materials_base"]; states=cfg["states"]
    mats_b=build_materials(base, states["bleached"]); mats_c=build_materials(base, states["colored"])
    Rb,Tb,Ab=tmm_RT(lam,layers,mats_b); Rc,Tc,Ac=tmm_RT(lam,layers,mats_c)
    Lt_b=luminous_transmittance(lam,Tb); Lt_c=luminous_transmittance(lam,Tc)
    np.savetxt(os.path.join(out,"window_spectra_bleached.csv"), np.column_stack([lam,Rb,Tb,Ab]), delimiter=",", header="lambda_um,R,T,A", comments="")
    np.savetxt(os.path.join(out,"window_spectra_colored.csv"),  np.column_stack([lam,Rc,Tc,Ac]), delimiter=",", header="lambda_um,R,T,A", comments="")
    prof=cfg["kinetics"]["voltage_profile"]; tau_on=cfg["kinetics"]["tau_s_color"]; tau_off=cfg["kinetics"]["tau_s_bleach"]
    t,V,X=simulate_kinetics(prof,tau_on,tau_off,dt=0.05); Lt_t=Lt_b*(1-X)+Lt_c*X
    np.savetxt(os.path.join(out,"window_kinetics.csv"), np.column_stack([t,V,X,Lt_t]), delimiter=",", header="t_s,Voltage_V,StateFractionColored,LuminousTransmittance", comments="")
    json.dump({"LuminousT_bleached":float(Lt_b),"LuminousT_colored":float(Lt_c),"Switch_time_est_color_s":tau_on,"Switch_time_est_bleach_s":tau_off},
              open(os.path.join(out,"window_summary.json"),"w"), indent=2)
    print("Window model done. Photopic T:", Lt_b, Lt_c)
if __name__=="__main__": main()
