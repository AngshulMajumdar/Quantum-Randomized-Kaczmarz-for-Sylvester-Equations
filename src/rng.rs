//! Stable 64-bit deterministic generator and Box–Muller Gaussian draws.
#[derive(Clone,Debug)]
pub struct Rng {state:u64}
impl Rng {
    pub fn new(seed:u64)->Self{Self{state:seed}}
    pub fn next_u64(&mut self)->u64{
        self.state=self.state.wrapping_add(0x9E3779B97F4A7C15);
        let mut z=self.state;
        z=(z^(z>>30)).wrapping_mul(0xBF58476D1CE4E5B9);
        z=(z^(z>>27)).wrapping_mul(0x94D049BB133111EB);
        z^(z>>31)
    }
    pub fn uniform(&mut self)->f64{(self.next_u64()>>11) as f64 * (1.0/(1u64<<53) as f64)}
    pub fn between(&mut self,a:f64,b:f64)->f64{a+(b-a)*self.uniform()}
    pub fn index(&mut self,n:usize)->usize{assert!(n>0);(self.next_u64() as u128*n as u128 >>64) as usize}
    pub fn normal(&mut self)->f64 {
        let u=1.0-self.uniform(); let v=self.uniform();
        (-2.0*u.ln()).sqrt()*(std::f64::consts::TAU*v).cos()
    }
}
