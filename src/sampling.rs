//! Exact norm-proportional sampling via two O(m+n) alias tables and rejection.
use crate::{matrix::Sylvester,rng::Rng};

#[derive(Clone,Debug)]
pub struct Alias {cutoff:Vec<f64>,alias:Vec<usize>}
impl Alias {
    pub fn new(weights:&[f64])->Result<Self,String>{
        let n=weights.len();
        if n==0 || weights.iter().any(|w|*w<0.0 || !w.is_finite()) {
            return Err("alias weights must be nonnegative and finite".into());
        }
        let total:f64=weights.iter().sum();
        if total<=0.0 {return Err("alias weights sum to zero".into());}
        let mut scaled=weights.iter().map(|w|w*n as f64/total).collect::<Vec<_>>();
        let mut small=Vec::new();let mut large=Vec::new();
        for (i,p) in scaled.iter().enumerate(){if *p<1.0{small.push(i);} else {large.push(i);}}
        let mut cutoff=vec![1.0;n];let mut alias=(0..n).collect::<Vec<_>>();
        while !small.is_empty() && !large.is_empty(){
            let s=small.pop().unwrap();let l=large.pop().unwrap();
            cutoff[s]=scaled[s];alias[s]=l;scaled[l]+=scaled[s]-1.0;
            if scaled[l]<1.0 {small.push(l);} else {large.push(l);}
        }
        Ok(Self{cutoff,alias})
    }
    pub fn sample(&self,rng:&mut Rng)->usize{
        let i=rng.index(self.cutoff.len());
        if rng.uniform()<self.cutoff[i] {i} else {self.alias[i]}
    }
}

#[derive(Clone,Debug)]
pub struct RowSampler {a_alias:Alias,b_alias:Alias,r2:Vec<f64>,s2:Vec<f64>,a_total:f64,b_total:f64}
impl RowSampler {
    pub fn new(eq:&Sylvester)->Result<Self,String>{
        let r2=(0..eq.m()).map(|i|(0..eq.m()).map(|k|eq.a.get(i,k).abs2()).sum()).collect::<Vec<f64>>();
        let s2=(0..eq.n()).map(|j|(0..eq.n()).map(|k|eq.b.get(k,j).abs2()).sum()).collect::<Vec<f64>>();
        let a_total=eq.n() as f64*r2.iter().sum::<f64>();
        let b_total=eq.m() as f64*s2.iter().sum::<f64>();
        let a_alias=Alias::new(&r2).or_else(|_|Alias::new(&vec![1.0;eq.m()]))?;
        let b_alias=Alias::new(&s2).or_else(|_|Alias::new(&vec![1.0;eq.n()]))?;
        if a_total+b_total<=0.0 {return Err("all coefficient rows are zero".into());}
        Ok(Self{a_alias,b_alias,r2,s2,a_total,b_total})
    }
    pub fn proposal(&self,eq:&Sylvester,rng:&mut Rng)->(usize,usize){
        if rng.uniform()*(self.a_total+self.b_total)<self.a_total {
            (self.a_alias.sample(rng),rng.index(eq.n()))
        } else { (rng.index(eq.m()),self.b_alias.sample(rng)) }
    }
    pub fn acceptance_rate(&self,eq:&Sylvester)->f64{eq.norm_sum_d()/(2.0*(self.a_total+self.b_total))}
    pub fn sample(&self,eq:&Sylvester,rng:&mut Rng)->Result<(usize,usize),String>{
        // The exact distribution conditional on acceptance is d_ij / D.
        // Check against invalid rows rather than silently biasing the distribution.
        for _ in 0..10_000_000 {
            let (i,j)=self.proposal(eq,rng);
            let denom=2.0*(self.r2[i]+self.s2[j]);
            if denom>0.0 {
                let acceptance=(eq.d(i,j)/denom).clamp(0.0,1.0);
                if rng.uniform()<acceptance{return Ok((i,j));}
            }
        }
        Err("sampler rejection budget exceeded: near-cancellation requires different access".into())
    }
    pub fn frozen(&self,eq:&Sylvester,rng:&mut Rng,t:usize)->Result<Vec<(usize,usize)>,String>{
        (0..t).map(|_|self.sample(eq,rng)).collect()
    }
}
