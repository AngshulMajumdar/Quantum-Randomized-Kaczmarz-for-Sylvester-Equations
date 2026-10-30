//! Row-major matrices and the Sylvester action. The solver does not form a Kronecker matrix.
use crate::complex::{C64, norm2};

#[derive(Clone, Debug)]
pub struct Matrix {
    pub rows: usize,
    pub cols: usize,
    pub data: Vec<C64>,
}
impl Matrix {
    pub fn zeros(rows: usize, cols: usize) -> Self {
        Self {rows,cols,data:vec![C64::ZERO;rows*cols]}
    }
    pub fn identity(n: usize) -> Self {
        let mut x=Self::zeros(n,n);
        for k in 0..n { x.set(k,k,C64::ONE); }
        x
    }
    pub fn get(&self,i:usize,j:usize)->C64 { self.data[i*self.cols+j] }
    pub fn set(&mut self,i:usize,j:usize,v:C64){ self.data[i*self.cols+j]=v; }
    pub fn add(&mut self,i:usize,j:usize,v:C64){ self.data[i*self.cols+j]+=v; }
    pub fn transpose(&self)->Self{
        let mut out=Self::zeros(self.cols,self.rows);
        for i in 0..self.rows { for j in 0..self.cols {out.set(j,i,self.get(i,j));} }
        out
    }
    pub fn adjoint(&self)->Self{
        let mut out=self.transpose();
        for z in &mut out.data { *z=z.conj(); }
        out
    }
    pub fn frob(&self)->f64 {norm2(&self.data).sqrt()}
    pub fn scaled(&self,scale:C64)->Self {Self {rows:self.rows,cols:self.cols,data:self.data.iter().map(|z|*z*scale).collect()}}
    pub fn mm(&self,b:&Self)->Self {
        assert_eq!(self.cols,b.rows);
        let mut out=Self::zeros(self.rows,b.cols);
        for i in 0..self.rows { for k in 0..self.cols {for j in 0..b.cols {
            out.add(i,j,self.get(i,k)*b.get(k,j));
        }}}
        out
    }
    pub fn diff_frob(&self,b:&Self)->f64 {
        assert_eq!((self.rows,self.cols),(b.rows,b.cols));
        self.data.iter().zip(&b.data).map(|(u,v)|(*u-*v).abs2()).sum::<f64>().sqrt()
    }
    pub fn is_zero(&self)->bool { self.data.iter().all(|z|*z==C64::ZERO) }
    pub fn vec_col(&self)->Vec<C64>{
        let mut out=Vec::with_capacity(self.rows*self.cols);
        for j in 0..self.cols {for i in 0..self.rows {out.push(self.get(i,j));}}
        out
    }
    pub fn from_col_vec(rows:usize,cols:usize,data:&[C64])->Self {
        assert_eq!(data.len(),rows*cols);
        let mut out=Self::zeros(rows,cols);
        for j in 0..cols {for i in 0..rows{out.set(i,j,data[j*rows+i]);}}
        out
    }
}

#[derive(Clone, Debug)]
pub struct Sylvester {pub a:Matrix,pub b:Matrix,pub c:Matrix}
impl Sylvester {
    pub fn new(a:Matrix,b:Matrix,c:Matrix)->Result<Self,String>{
        if a.rows==0||b.rows==0||a.rows!=a.cols||b.rows!=b.cols
            || c.rows!=a.rows||c.cols!=b.rows {
            return Err("expected square nonempty A, B, with C of size m x n".into());
        }
        Ok(Self{a,b,c})
    }
    pub fn m(&self)->usize{self.a.rows}
    pub fn n(&self)->usize{self.b.rows}
    pub fn apply(&self,x:&Matrix)->Matrix {
        let mut ax=self.a.mm(x);
        let xb=x.mm(&self.b);
        for (u,v) in ax.data.iter_mut().zip(&xb.data) {*u+=*v;}
        ax
    }
    pub fn residual(&self,x:&Matrix)->f64 {self.apply(x).diff_frob(&self.c)}
    pub fn forcing_from(a:Matrix,b:Matrix,x:&Matrix)->Self{
        let zeros=Matrix::zeros(a.rows,b.rows);
        let mut eq=Self::new(a,b,zeros).expect("valid dimensions");
        eq.c=eq.apply(x);
        eq
    }
    /// A†e_i, as a vector; entries are conjugated row i of A.
    pub fn a_row_conj(&self,i:usize)->Vec<C64>{(0..self.m()).map(|k|self.a.get(i,k).conj()).collect()}
    /// b_j=Be_j, as a vector.
    pub fn b_column(&self,j:usize)->Vec<C64>{(0..self.n()).map(|k|self.b.get(k,j)).collect()}
    pub fn d(&self,i:usize,j:usize)->f64{
        let r2=(0..self.m()).map(|k|self.a.get(i,k).abs2()).sum::<f64>();
        let s2=(0..self.n()).map(|k|self.b.get(k,j).abs2()).sum::<f64>();
        let cross=2.0*(self.a.get(i,i)*self.b.get(j,j).conj()).re;
        r2+s2+cross
    }
    /// Normalized vec(Γ_ij), with column stacking and the required complex conjugates.
    pub fn normalized_row(&self,i:usize,j:usize)->Result<Vec<C64>,String>{
        let d=self.d(i,j);
        if d<=0.0 || !d.is_finite(){return Err(format!("row ({i},{j}) is zero or invalid"));}
        let mut g=vec![C64::ZERO;self.m()*self.n()];
        for k in 0..self.m(){g[j*self.m()+k]+=self.a.get(i,k).conj();}
        for l in 0..self.n(){g[l*self.m()+i]+=self.b.get(l,j).conj();}
        for z in &mut g { *z=*z/d.sqrt(); }
        Ok(g)
    }
    pub fn beta(&self,i:usize,j:usize)->C64{self.c.get(i,j)/self.d(i,j).sqrt()}
    pub fn norm_sum_d(&self)->f64{
        let mut total=0.0;
        for i in 0..self.m(){for j in 0..self.n(){total+=self.d(i,j);}}
        total
    }
    /// Dense Kronecker direct solver is used only for small validation systems, never in QRK.
    pub fn direct_solve(&self)->Result<Matrix,String>{
        let m=self.m();let n=self.n();let size=m*n;
        if size>4096 {return Err("dense reference solve restricted to mn <= 4096".into());}
        let mut aug=vec![vec![C64::ZERO;size+1];size];
        for j in 0..n {for i in 0..m {
            let r=j*m+i;
            for k in 0..m {aug[r][j*m+k]+=self.a.get(i,k);}
            for l in 0..n {aug[r][l*m+i]+=self.b.get(l,j);}
            aug[r][size]=self.c.get(i,j);
        }}
        for p in 0..size {
            let mut best=p;
            for r in p+1..size {if aug[r][p].abs2()>aug[best][p].abs2(){best=r;}}
            if aug[best][p].abs()<1e-13 {return Err("Sylvester operator is singular or numerically unresolved".into());}
            aug.swap(p,best);
            let pivot=aug[p][p];
            for c in p..=size {aug[p][c]=aug[p][c]/pivot;}
            for r in 0..size {if r!=p {
                let factor=aug[r][p];
                if factor!=C64::ZERO {
                    for c in p..=size {let v=aug[p][c]; aug[r][c]-=factor*v;}
                }
            }}
        }
        let sol=(0..size).map(|i|aug[i][size]).collect::<Vec<_>>();
        Ok(Matrix::from_col_vec(m,n,&sol))
    }
}
