//! QRK, matrix-sign, spectral, Lyapunov-integral and mean-channel numerical kernels.
//! Baseline kernels represent classical emulations of published mathematical maps,
//! NOT implementations of the authors' quantum block encodings or fault-tolerant circuits.
use crate::complex::{C64, inner};
use crate::matrix::{Matrix,Sylvester};
use crate::rng::Rng;
use crate::sampling::RowSampler;
use crate::linalg::{add_scaled, scale, inverse, expm, real_symmetric_eigen,
    hermitian_part_positive,is_normal,is_positive_semidefinite};

#[derive(Clone,Debug)]
pub struct QrkOutput {
    pub x:Matrix, pub frozen_sequence:Vec<(usize,usize)>,
    pub active_probability:f64, pub reference_probability:f64,
    pub normalization_log:f64, pub sigma:f64,
}

/// Exact paper Kaczmarz iterates with one frozen sample sequence and analytical
/// active-clock amplitudes. R is a supplied bound on ||X*||_F.
pub fn qrk(eq:&Sylvester,steps:usize,bound:f64,seed:u64)->Result<QrkOutput,String>{
    if eq.c.is_zero(){return Ok(QrkOutput {x:Matrix::zeros(eq.m(),eq.n()),frozen_sequence:vec![],active_probability:1.0,reference_probability:1.0,normalization_log:0.0,sigma:0.0});}
    if steps==0||!bound.is_finite()||bound<=0.0{return Err("QRK requires steps>0 and a positive R bound".into());}
    let sigma=(steps as f64).sqrt()*bound;
    let sampler=RowSampler::new(eq)?;
    let mut rng=Rng::new(seed);
    let tau=sampler.frozen(eq,&mut rng,steps)?;
    let mut x=vec![C64::ZERO;eq.m()*eq.n()];
    let mut log_alpha=0.0;
    for &(i,j) in &tau {
        let g=eq.normalized_row(i,j)?;
        let beta=eq.beta(i,j);
        let error=beta-inner(&g,&x);
        for (u,v) in x.iter_mut().zip(g.iter()){*u+=*v*error;}
        log_alpha+=0.5*(1.0+(beta/sigma).abs2()).ln();
    }
    let xmat=Matrix::from_col_vec(eq.m(),eq.n(),&x);
    let active=(1.0+(xmat.frob()/sigma).powi(2))*(-2.0*log_alpha).exp();
    let reference=(-2.0*log_alpha).exp();
    Ok(QrkOutput{x:xmat,frozen_sequence:tau,active_probability:active,reference_probability:reference,normalization_log:log_alpha,sigma})
}
/// Generic vectorized norm-proportional Kaczmarz with explicit K in memory.
/// This is a classical ideal-row-oracle comparator (same projections as QRK).
pub fn generic_kaczmarz(eq:&Sylvester,sequence:&[(usize,usize)])->Result<Matrix,String>{
    let m=eq.m();let n=eq.n();let sz=m*n;
    if sz>1024{return Err("explicit Kronecker baseline restricted to mn <= 1024".into());}
    let mut rows=vec![vec![C64::ZERO;sz];sz];
    for j in 0..n {for i in 0..m {
        let idx=j*m+i;
        for k in 0..m {rows[idx][j*m+k]+=eq.a.get(i,k);}
        for l in 0..n {rows[idx][l*m+i]+=eq.b.get(l,j);}
    }}
    let mut x=vec![C64::ZERO;sz];
    for &(i,j) in sequence {
        let row=&rows[j*m+i];
        let norm=row.iter().map(|z|z.abs2()).sum::<f64>();
        if norm<=0.0 {return Err("zero Kronecker row".into());}
        let predicted=row.iter().zip(x.iter()).map(|(a,b)|*a*(*b)).fold(C64::ZERO,|a,b|a+b);
        let update=(eq.c.get(i,j)-predicted)/norm;
        for (u,v) in x.iter_mut().zip(row.iter()){*u+=v.conj()*update;}
    }
    Ok(Matrix::from_col_vec(m,n,&x))
}
/// Somma–Low–Berry–Babbush mathematical solution-map emulator.
/// (i) real-symmetric spectral functional calculus, including indefinite spectra;
/// (ii) Laplace functional calculus on strongly positive-real A,B.
/// Does NOT simulate the paper's block-encoding, QSVT or query complexity.
pub fn somma(eq:&Sylvester)->Result<Matrix,String>{
    // Complex diagonal normal case: diagonal functional calculus needs no
    // real-symmetric restriction and retains complex phases exactly.
    let diagonal = |a:&Matrix| -> bool {
        (0..a.rows).all(|i|(0..a.cols).all(|j|i==j || a.get(i,j).abs()<1e-13))
    };
    if diagonal(&eq.a)&&diagonal(&eq.b){
        let mut x=Matrix::zeros(eq.m(),eq.n());
        for i in 0..eq.m(){for j in 0..eq.n(){
            let denom=eq.a.get(i,i)+eq.b.get(j,j);
            if denom.abs()<1e-13{return Err("Somma diagonal map: intersecting spectra".into());}
            x.set(i,j,eq.c.get(i,j)/denom);
        }}
        return Ok(x)
    }
    if let (Ok((da,ua)),Ok((db,ub)))=(real_symmetric_eigen(&eq.a),real_symmetric_eigen(&eq.b)) {
        let transformed=ua.adjoint().mm(&eq.c).mm(&ub);
        let mut x=Matrix::zeros(eq.m(),eq.n());
        for i in 0..eq.m(){for j in 0..eq.n(){
            let denom=da[i]+db[j];
            if denom.abs()<1e-11{return Err("Somma spectral map: intersecting Sylvester spectra".into());}
            x.set(i,j,transformed.get(i,j)/denom);
        }}
        return Ok(ua.mm(&x).mm(&ub.adjoint()));
    }
    if !(hermitian_part_positive(&eq.a,1.0)&&hermitian_part_positive(&eq.b,1.0)) {
        return Err("Somma emulator: no supported symmetric or positive-real spectral regime".into());
    }
    // X=integral_0^infty exp(-At) C exp(-Bt) dt.
    // Composite trapezoid with explicit step and tail checks, independent of direct solve.
    let h=0.02;
    let ea=expm(&scale(&eq.a,-h));let eb=expm(&scale(&eq.b,-h));
    let mut q=eq.c.clone();let mut sum=scale(&q,0.5*h);
    for k in 1..=15000 {
        q=ea.mm(&q).mm(&eb);
        let weight=if k==15000 {0.5} else {1.0};
        sum=add_scaled(&sum,&q,C64::from(weight*h));
        if k>10 && q.frob()*h<1e-12*eq.c.frob().max(1e-14){break;}
    }
    Ok(sum)
}

/// Wang–Liu sign-embedding matrix-function emulator.
/// For separating half planes sign([[A,C],[0,-B]])_12 = 2s X.
/// Newton matrix-sign iterations are classical matrix operations, not quantum
/// log-sinc approximation, scaled multiplexing, or quantum shifted inverses.
pub fn wang_liu(eq:&Sylvester)->Result<Matrix,String>{
    let s=if hermitian_part_positive(&eq.a,1.0)&&hermitian_part_positive(&eq.b,1.0){1.0}
        else if hermitian_part_positive(&eq.a,-1.0)&&hermitian_part_positive(&eq.b,-1.0){-1.0}
        else {return Err("Wang–Liu: no certified common half-plane separation".into());};
    let m=eq.m();let n=eq.n();let mut z=Matrix::zeros(m+n,m+n);
    for i in 0..m {for j in 0..m{z.set(i,j,eq.a.get(i,j));}}
    for i in 0..n {for j in 0..n{z.set(m+i,m+j,-eq.b.get(i,j));}}
    for i in 0..m {for j in 0..n {z.set(i,m+j,eq.c.get(i,j));}}
    let mut converged=false;
    for _ in 0..120 {
        let inv=inverse(&z)?;
        let next=scale(&add_scaled(&z,&inv,C64::ONE),0.5);
        if next.diff_frob(&z)<1e-11*(1.0+next.frob()) {z=next;converged=true;break;}
        z=next;
    }
    if !converged{return Err("matrix sign iteration did not converge".into());}
    let mut x=Matrix::zeros(m,n);
    for i in 0..m {for j in 0..n{x.set(i,j,z.get(i,m+j)/(2.0*s));}}
    Ok(x)
}

fn lyapunov_data(eq:&Sylvester,require_normal:bool)->Result<Matrix,String>{
    if eq.m()!=eq.n()||eq.b.diff_frob(&eq.a.adjoint())>1e-9*(1.0+eq.a.frob()) {
        return Err("Lyapunov baseline requires B=A*".into());
    }
    if !hermitian_part_positive(&eq.a,-1.0){return Err("Lyapunov baseline requires certified Hurwitz A".into());}
    if require_normal && !is_normal(&eq.a,1e-9){return Err("Benedetti baseline requires normal A".into());}
    let omega=scale(&eq.c,-1.0);
    if !is_positive_semidefinite(&omega,1e-9){return Err("Lyapunov baseline requires positive-semidefinite -C".into());}
    Ok(omega)
}
/// Clayton et al. Lyapunov subroutine's deterministic integral-value emulator.
/// Does not implement quantum policy gradient or its quantum access primitives.
pub fn clayton(eq:&Sylvester)->Result<Matrix,String>{
    let omega=lyapunov_data(eq,false)?;
    let dt=0.02;let evolution=expm(&scale(&eq.a,dt));
    let mut q=omega;let mut x=scale(&q,0.5*dt);
    for k in 1..=2000 {
        q=evolution.mm(&q).mm(&evolution.adjoint());
        x=add_scaled(&x,&q,C64::from(dt));
        if k>10 && q.frob()*dt<1e-12*x.frob().max(1e-14){break;}
    }
    Ok(x)
}
/// Benedetti–Rosmanis–Rosenkranz expected mixed-state channel emulator.
/// E=e^{A dt}; E_{k+1}=E E_k E† with positive semidefinite source.
/// Reports the unnormalized expectation X (not an actual stochastic quantum circuit).
pub fn benedetti(eq:&Sylvester)->Result<Matrix,String>{
    let omega=lyapunov_data(eq,true)?;
    let dt=0.01;
    let e=expm(&scale(&eq.a,dt));
    let mut q=omega;let mut x=scale(&q,0.5*dt);
    for k in 1..=4000 {
        q=e.mm(&q).mm(&e.adjoint());
        x=add_scaled(&x,&q,C64::from(dt));
        if k>10 && dt*q.frob()<1e-12*x.frob().max(1e-14){break;}
    }
    Ok(x)
}

/// Entrywise phase-sensitive finite-shot measurement of a normalized augmented
/// state. Each Re/Im coherence is inferred from three-outcome +/-/other POVMs.
/// Per-outcome depolarization and symmetric readout attenuation are mitigated.
pub fn finite_shots(x:&Matrix,sigma:f64,shots:usize,seed:u64,visibility:f64,bitflip:f64)->Result<Matrix,String>{
    if sigma<=0.0||shots==0||visibility<=0.0||visibility>1.0||bitflip<0.0||bitflip>=0.5 {
        return Err("invalid finite-shot parameters".into());
    }
    let mut rng=Rng::new(seed);
    let norm=x.frob();let z=sigma*sigma+norm*norm;
    let p_ref=sigma*sigma/z;
    let attenuation=visibility*(1.0-2.0*bitflip);
    let mut out=Matrix::zeros(x.rows,x.cols);
    for k in 0..x.data.len() {
        let amp=x.data[k]*C64::from(sigma/z);
        let density=(x.data[k].abs2()+sigma*sigma)/z;
        let estimate = |signal:f64, rng:&mut Rng| -> f64 {
            let plus=visibility*0.5*(density+signal)+(1.0-visibility)/((x.data.len()+1) as f64);
            // Visibility noise is treated as uniform depolarization on the augmented space.
            let minus=visibility*0.5*(density-signal)+(1.0-visibility)/((x.data.len()+1) as f64);
            let plus=plus.clamp(0.0,1.0);let minus=minus.clamp(0.0,1.0-plus);
            let mut acc=0.0;
            for _ in 0..shots {
                let u=rng.uniform();
                let outcome=if u<plus {1.0} else if u<plus+minus {-1.0} else {0.0};
                acc+=if outcome!=0.0 && rng.uniform()<bitflip{-outcome}else{outcome};
            }
            acc/(shots as f64)/attenuation
        };
        let re=estimate(2.0*amp.re,&mut rng)/2.0;
        let im=estimate(2.0*amp.im,&mut rng)/2.0;
        out.data[k]=C64::new(re,im)*C64::from(sigma/p_ref);
    }
    Ok(out)
}
