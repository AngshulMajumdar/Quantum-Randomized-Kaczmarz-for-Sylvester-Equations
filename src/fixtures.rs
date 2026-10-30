//! Four paper instance families and two additional diagnostic families.
use crate::complex::C64;
use crate::matrix::{Matrix,Sylvester};
use crate::rng::Rng;

#[derive(Clone,Copy,Debug,PartialEq,Eq)]
pub enum Family {Screened,Convection,Wave,Control,NormalLyapunov,Complex}
impl Family {
    pub fn all()->[Self;6]{[Self::Screened,Self::Convection,Self::Wave,Self::Control,Self::NormalLyapunov,Self::Complex]}
    pub fn paper()->[Self;4]{[Self::Screened,Self::Convection,Self::Wave,Self::Control]}
    pub fn as_str(self)->&'static str{match self {
        Self::Screened=>"screened",Self::Convection=>"convection",Self::Wave=>"wave",
        Self::Control=>"control",Self::NormalLyapunov=>"normal-lyapunov",Self::Complex=>"complex"
    }}
    pub fn parse(s:&str)->Option<Self>{Self::all().into_iter().find(|f|f.as_str()==s)}
}
fn tridiag(n:usize,lo:f64,diag:f64,hi:f64)->Matrix{
    let mut a=Matrix::zeros(n,n);
    for i in 0..n {
        a.set(i,i,C64::from(diag));
        if i>0 {a.set(i,i-1,C64::from(lo));}
        if i+1<n{a.set(i,i+1,C64::from(hi));}
    }
    a
}
fn random_x(n:usize,rng:&mut Rng,complex:bool)->Matrix {
    let mut x=Matrix::zeros(n,n);
    for i in 0..n {for j in 0..n {
        x.set(i,j,C64::new(rng.normal(),if complex{rng.normal()}else{0.0}));
    }}x
}
pub fn make(f:Family,n:usize,seed:u64)->Result<Sylvester,String>{
    if n<2 {return Err("fixture requires n >= 2".into());}
    if f==Family::Wave && n%2!=0{return Err("wave family requires even n".into());}
    let mut r=Rng::new(seed);
    let (a,b,c_opt)=match f {
        Family::Screened=>{
            (tridiag(n,-1.0,2.4+r.between(0.0,0.25),-1.0),
             tridiag(n,-0.75,1.9+r.between(0.0,0.25),-0.75),None)
        },
        Family::Convection=>{
            // Draw parameters once for each coefficient, as in the paper
            let ma=r.between(0.18,0.35);let mb=r.between(0.15,0.30);
            let va=r.between(0.5,1.2);let vb=r.between(0.4,1.0);
            let a=tridiag(n,-ma-va,2.0*ma+va+0.35,-ma);
            let b=tridiag(n,-mb-vb,2.0*mb+vb+0.25,-mb);
            (a,b,None)
        },
        Family::Wave=>{
            let mut a=Matrix::zeros(n,n);let mut b=Matrix::zeros(n,n);
            for k in (0..n).step_by(2){
                let va=r.between(2.0,2.4);let vb=r.between(0.9,1.3);
                a.set(k,k+1,C64::from(va));a.set(k+1,k,C64::from(va));
                b.set(k,k+1,C64::from(vb));b.set(k+1,k,C64::from(vb));
            }(a,b,None)
        },
        Family::Control=>{
            let gamma=r.between(0.45,0.75);
            let mut a=Matrix::identity(n).scaled(C64::from(-1.0));
            for i in 0..n-1{a.set(i,i+1,C64::from(gamma));}
            let b=a.transpose();
            let mut g=Matrix::zeros(n,n);
            for z in &mut g.data {*z=C64::from(r.normal());}
            let mut omega=g.mm(&g.transpose()).scaled(C64::from(1.0/n as f64));
            for i in 0..n {omega.add(i,i,C64::from(0.25));}
            (a,b,Some(omega.scaled(C64::from(-1.0))))
        },
        Family::NormalLyapunov=>{
            let mut a=Matrix::zeros(n,n);
            for i in 0..n {a.set(i,i,C64::new(-1.0-r.uniform(),0.2*(i as f64)));}
            let b=a.adjoint();
            let mut g=Matrix::zeros(n,n);
            for z in &mut g.data{*z=C64::new(r.normal(),r.normal());}
            let mut omega=g.mm(&g.adjoint()).scaled(C64::from(1.0/(n as f64)));
            for i in 0..n{omega.add(i,i,C64::from(0.25));}
            (a,b,Some(omega.scaled(C64::from(-1.0))))
        },
        Family::Complex=>{
            let mut a=Matrix::zeros(n,n);let mut b=Matrix::zeros(n,n);
            for i in 0..n {
                a.set(i,i,C64::new(3.0+r.uniform(),r.between(-0.2,0.2)));
                b.set(i,i,C64::new(2.0+r.uniform(),r.between(-0.2,0.2)));
                if i+1<n {
                    a.set(i,i+1,C64::new(-0.8,r.between(-0.3,0.3)));
                    a.set(i+1,i,C64::new(-0.5,r.between(-0.3,0.3)));
                    b.set(i,i+1,C64::new(-0.6,r.between(-0.3,0.3)));
                    b.set(i+1,i,C64::new(-0.4,r.between(-0.3,0.3)));
                }
            }(a,b,None)
        }
    };
    if let Some(c)=c_opt {Sylvester::new(a,b,c)} else {
        let x=random_x(n,&mut r,f==Family::Complex);
        Ok(Sylvester::forcing_from(a,b,&x))
    }
}
