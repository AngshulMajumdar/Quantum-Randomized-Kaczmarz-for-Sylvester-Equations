//! Dependency-free double-precision complex arithmetic.
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub, SubAssign};

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}
impl C64 {
    pub const ZERO: Self = Self { re: 0.0, im: 0.0 };
    pub const ONE: Self = Self { re: 1.0, im: 0.0 };
    pub fn new(re: f64, im: f64) -> Self { Self { re, im } }
    pub fn conj(self) -> Self { Self::new(self.re, -self.im) }
    pub fn abs2(self) -> f64 { self.re * self.re + self.im * self.im }
    pub fn abs(self) -> f64 { self.abs2().sqrt() }
    pub fn scale(self, s: f64) -> Self { Self::new(self.re * s, self.im * s) }
}
impl From<f64> for C64 { fn from(x: f64) -> Self { Self::new(x, 0.0) } }
impl Add for C64 { type Output=Self; fn add(self, x: Self) -> Self { Self::new(self.re+x.re, self.im+x.im) } }
impl Sub for C64 { type Output=Self; fn sub(self, x: Self) -> Self { Self::new(self.re-x.re, self.im-x.im) } }
impl Neg for C64 { type Output=Self; fn neg(self) -> Self { Self::new(-self.re,-self.im) } }
impl Mul for C64 {
    type Output=Self;
    fn mul(self, x: Self) -> Self {
        Self::new(self.re*x.re-self.im*x.im,self.re*x.im+self.im*x.re)
    }
}
impl Div for C64 {
    type Output=Self;
    fn div(self, x: Self) -> Self {
        let d=x.abs2();
        Self::new((self.re*x.re+self.im*x.im)/d,(self.im*x.re-self.re*x.im)/d)
    }
}
impl Mul<f64> for C64 { type Output=Self; fn mul(self, x: f64)->Self { self.scale(x) } }
impl Div<f64> for C64 { type Output=Self; fn div(self, x: f64)->Self { self.scale(1.0/x) } }
impl AddAssign for C64 { fn add_assign(&mut self, x: Self) { *self=*self+x; } }
impl SubAssign for C64 { fn sub_assign(&mut self, x: Self) { *self=*self-x; } }

pub fn inner(x: &[C64], y: &[C64]) -> C64 {
    assert_eq!(x.len(),y.len());
    x.iter().zip(y).fold(C64::ZERO,|acc,(a,b)|acc+a.conj()*(*b))
}
pub fn norm2(x: &[C64]) -> f64 { x.iter().map(|v|v.abs2()).sum() }
