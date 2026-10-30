use qrk_sylv::complex::C64;
use qrk_sylv::fixtures::{make,Family};
use qrk_sylv::matrix::{Matrix,Sylvester};
use qrk_sylv::solvers::{qrk,generic_kaczmarz,somma,wang_liu,clayton,benedetti,finite_shots};
use qrk_sylv::clock::simulate;
use qrk_sylv::linalg::{inverse,expm,real_symmetric_eigen};

fn relative(a:&Matrix,b:&Matrix)->f64{a.diff_frob(b)/b.frob().max(1e-15)}

#[test]
fn complex_row_norm_and_projection(){
    let mut a=Matrix::identity(2);let mut b=Matrix::identity(2);
    a.set(0,0,C64::new(2.0,1.0));b.set(0,0,C64::new(3.0,2.0));
    let eq=Sylvester::new(a,b,Matrix::identity(2)).unwrap();
    // g_00 diag entry = conj(A_00+B_00) = 5 - 3i; d_00 = 34.
    assert!((eq.d(0,0)-34.0).abs()<1e-12);
    let row=eq.normalized_row(0,0).unwrap();
    let n=row.iter().map(|z|z.abs2()).sum::<f64>();
    assert!((n-1.0).abs()<1e-12);
}
#[test]
fn linear_algebra_basics(){
    let a=make(Family::Complex,4,918).unwrap().a;
    let inv=inverse(&a).unwrap();
    assert!(relative(&a.mm(&inv),&Matrix::identity(4))<1e-10);
    let ex=expm(&a.scaled(C64::from(0.01)));
    let em=expm(&a.scaled(C64::from(-0.01)));
    assert!(relative(&ex.mm(&em),&Matrix::identity(4))<1e-10);
    let symmetric=make(Family::Wave,4,918).unwrap().a;
    let (values,vecs)=real_symmetric_eigen(&symmetric).unwrap();
    let mut diagonal=Matrix::zeros(4,4);
    for i in 0..4{diagonal.set(i,i,C64::from(values[i]));}
    assert!(relative(&vecs.mm(&diagonal).mm(&vecs.adjoint()),&symmetric)<1e-10);
}
#[test]
fn qrk_matches_explicit_kronecker(){
    for family in [Family::Screened,Family::Convection,Family::Complex,Family::Wave] {
        let eq=make(family,2,123).unwrap();
        let truth=eq.direct_solve().unwrap();
        let output=qrk(&eq,400,truth.frob(),33).unwrap();
        let generic=generic_kaczmarz(&eq,&output.frozen_sequence).unwrap();
        assert!(relative(&output.x,&generic)<1e-11,"{}",family.as_str());
        assert!(output.active_probability>=(-1.0f64).exp()-1e-12);
    }
}
#[test]
fn actual_failure_clock_is_unitary_and_matches_qrk(){
    for f in [Family::Screened,Family::Complex] {
        let eq=make(f,2,111).unwrap();let truth=eq.direct_solve().unwrap();
        let qr=qrk(&eq,100,truth.frob(),423).unwrap();
        let clock=simulate(&eq,&qr.frozen_sequence,truth.frob()).unwrap();
        assert!(relative(&clock.reconstructed,&qr.x)<1e-11);
        assert!((clock.state_norm-1.0).abs()<1e-10);
        assert!((clock.active_probability-qr.active_probability).abs()<1e-10);
    }
}
#[test]
fn somma_spectral_and_laplace_maps(){
    for f in [Family::Screened,Family::Wave,Family::Convection,Family::NormalLyapunov] {
        let eq=make(f,4,37).unwrap();let truth=eq.direct_solve().unwrap();
        assert!(relative(&somma(&eq).unwrap(),&truth)<1e-3,"{}",f.as_str());
    }
}
#[test]
fn wang_sign_positive_and_negative_halfplanes(){
    for f in [Family::Screened,Family::Convection,Family::Control] {
        let eq=make(f,4,35).unwrap();let truth=eq.direct_solve().unwrap();
        assert!(relative(&wang_liu(&eq).unwrap(),&truth)<1e-8,"{}",f.as_str());
    }
    assert!(wang_liu(&make(Family::Wave,4,42).unwrap()).is_err());
}
#[test]
fn lyapunov_implementations_and_domain_checks(){
    let control=make(Family::Control,4,9).unwrap();
    let truth=control.direct_solve().unwrap();
    assert!(relative(&clayton(&control).unwrap(),&truth)<1e-3);
    assert!(benedetti(&control).is_err());
    let normal=make(Family::NormalLyapunov,4,9).unwrap();
    let truth=normal.direct_solve().unwrap();
    assert!(relative(&benedetti(&normal).unwrap(),&truth)<1e-3);
    assert!(relative(&clayton(&normal).unwrap(),&truth)<1e-3);
    assert!(clayton(&make(Family::Wave,4,9).unwrap()).is_err());
}
#[test]
fn finite_shots_exposes_complex_coherences(){
    let mut x=Matrix::zeros(2,2);x.set(0,1,C64::new(0.7,-0.8));
    let approx=finite_shots(&x,10.0,150_000,123,0.975,0.015).unwrap();
    assert!((approx.get(0,1).re-0.7).abs()<0.35);
    assert!((approx.get(0,1).im+0.8).abs()<0.35);
}

#[test]
fn alias_rejection_matches_exact_row_probabilities(){
    use qrk_sylv::rng::Rng;
    use qrk_sylv::sampling::RowSampler;
    let eq=make(Family::Complex,2,456).unwrap();
    let sampler=RowSampler::new(&eq).unwrap();
    let mut rng=Rng::new(917);
    let mut counts=[0usize;4];
    let trials=60_000;
    for _ in 0..trials {
        let (i,j)=sampler.sample(&eq,&mut rng).unwrap();
        counts[2*j+i]+=1;
    }
    let total=eq.norm_sum_d();
    for j in 0..2 {for i in 0..2 {
        let observed=counts[j*2+i] as f64/trials as f64;
        let exact=eq.d(i,j)/total;
        assert!((observed-exact).abs()<0.015, "({i},{j}): {observed} vs {exact}");
    }}
}
