//! Exact Julia-dilation state-vector emulator with distinct failure-clock sectors.
//! This allocates O((T+1)(mn+1)) complex amplitudes; it is for small validation.
use crate::complex::{C64,inner,norm2};
use crate::matrix::{Matrix,Sylvester};

#[derive(Clone,Debug)]
pub struct ClockOutput {
    pub active:Vec<C64>, pub failure:Vec<Vec<C64>>,
    pub active_probability:f64, pub state_norm:f64, pub reconstructed:Matrix,
}

/// Apply the 2-sector Julia dilation to the active sector at each step.
/// The failure sector for a step starts in zero and is not revisited.
pub fn simulate(eq:&Sylvester,sequence:&[(usize,usize)],bound:f64)->Result<ClockOutput,String>{
    let sz=eq.m()*eq.n();let t=sequence.len();
    if t==0||bound<=0.0||!bound.is_finite(){return Err("clock simulation requires nonzero steps and positive bound".into());}
    let sigma=(t as f64).sqrt()*bound;
    let mut active=vec![C64::ZERO;sz+1];
    active[sz]=C64::ONE;
    let mut failure=Vec::with_capacity(t);
    let mut cumulative=1.0;
    for &(i,j) in sequence {
        let g=eq.normalized_row(i,j)?;
        let beta=eq.beta(i,j);
        let b=beta/sigma;
        let alpha=(1.0+b.abs2()).sqrt();
        let proj=inner(&g,&active[0..sz]);
        let delta=b.abs()/alpha;
        let mut rejected=vec![C64::ZERO;sz+1];
        for k in 0..sz {
            rejected[k]=g[k]*proj + (active[k]-g[k]*proj)*delta;
        }
        // D=(I-K†K)^(1/2) on the source active state.
        // The operator is P_g + delta (I_data - P_g), with zero reference part.
        let mut next=vec![C64::ZERO;sz+1];
        for k in 0..sz{
            next[k]=(active[k]-g[k]*proj+g[k]*b*active[sz])/alpha;
        }
        next[sz]=active[sz]/alpha;
        cumulative*=alpha;
        active=next;
        failure.push(rejected);
    }
    let p=norm2(&active);
    let norm=p+failure.iter().map(|v|norm2(v)).sum::<f64>();
    let x=Matrix::from_col_vec(eq.m(),eq.n(),&active[0..sz].iter().map(|z|(*z)*sigma*cumulative).collect::<Vec<_>>());
    Ok(ClockOutput{active,failure,active_probability:p,state_norm:norm,reconstructed:x})
}
