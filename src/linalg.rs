//! Small dense complex linear algebra for reference algorithms (no external crates).
use crate::complex::C64;
use crate::matrix::Matrix;

pub fn add_scaled(a:&Matrix,b:&Matrix,s:C64)->Matrix {
    assert_eq!((a.rows,a.cols),(b.rows,b.cols));
    let mut x=a.clone();
    for (u,v) in x.data.iter_mut().zip(b.data.iter()) {*u+=*v*s;}
    x
}
pub fn scale(a:&Matrix,s:f64)->Matrix {a.scaled(C64::from(s))}
pub fn inverse(a:&Matrix)->Result<Matrix,String>{
    if a.rows!=a.cols {return Err("inverse requires square matrix".into());}
    let n=a.rows;
    let mut x=a.clone();let mut inv=Matrix::identity(n);
    for k in 0..n {
        let mut p=k;
        for i in k+1..n {if x.get(i,k).abs2()>x.get(p,k).abs2(){p=i;}}
        if x.get(p,k).abs()<1e-14 {return Err("matrix singular in inverse".into());}
        for j in 0..n {
            let t=x.get(k,j); x.set(k,j,x.get(p,j));x.set(p,j,t);
            let t=inv.get(k,j); inv.set(k,j,inv.get(p,j));inv.set(p,j,t);
        }
        let d=x.get(k,k);
        for j in 0..n {x.set(k,j,x.get(k,j)/d);inv.set(k,j,inv.get(k,j)/d);}
        for i in 0..n {if i!=k {
            let f=x.get(i,k);
            if f.abs2()>0.0 {for j in 0..n {
                x.set(i,j,x.get(i,j)-f*x.get(k,j));
                inv.set(i,j,inv.get(i,j)-f*inv.get(k,j));
            }}
        }}
    }
    Ok(inv)
}
/// Scaling-and-squaring Taylor action for small dense matrices.
pub fn expm(a:&Matrix)->Matrix {
    assert_eq!(a.rows,a.cols);
    let n=a.rows;
    let mut scale_pow=0usize;
    let mut scaled=a.clone();
    while scaled.frob()>0.5 && scale_pow<40 {scaled=scale(&scaled,0.5);scale_pow+=1;}
    let mut result=Matrix::identity(n);
    let mut term=Matrix::identity(n);
    for k in 1..=80 {
        term=scale(&term.mm(&scaled),1.0/(k as f64));
        result=add_scaled(&result,&term,C64::ONE);
        if term.frob()<1e-15 {break;}
    }
    for _ in 0..scale_pow {result=result.mm(&result);}
    result
}
pub fn hermitian_part_positive(a:&Matrix,sign:f64)->bool {
    if a.rows!=a.cols{return false;}
    // Strict diagonal dominance of the signed Hermitian part is a sufficient,
    // deliberately conservative half-plane certificate.
    for i in 0..a.rows {
        let d=sign*a.get(i,i).re;
        let mut sum=0.0;
        for j in 0..a.cols {if i!=j {
            let h=(a.get(i,j)+a.get(j,i).conj())*0.5;
            sum+=h.abs();
        }}
        if d<=sum+1e-10 {return false;}
    }
    true
}
pub fn is_hermitian(a:&Matrix,tol:f64)->bool {
    a.rows==a.cols && a.diff_frob(&a.adjoint())<=tol*(1.0+a.frob())
}
pub fn is_normal(a:&Matrix,tol:f64)->bool{
    let d=a.mm(&a.adjoint()).diff_frob(&a.adjoint().mm(a));
    d<=tol*(1.0+a.frob()*a.frob())
}
pub fn is_positive_semidefinite(a:&Matrix,tol:f64)->bool {
    if !is_hermitian(a,tol){return false;}
    let n=a.rows;
    // Complex LDL* without pivoting; shift by tolerance for semidefinite case.
    let mut l=Matrix::identity(n);let mut d=vec![0.0f64;n];
    for i in 0..n {
        let mut s=a.get(i,i).re+tol;
        for k in 0..i {s-=l.get(i,k).abs2()*d[k];}
        if s < -tol*10.0{return false;}
        d[i]=s.max(tol);
        for j in i+1..n {
            let mut v=a.get(j,i);
            for k in 0..i {v-=l.get(j,k)*l.get(i,k).conj()*d[k];}
            l.set(j,i,v/d[i]);
        }
    }
    true
}
/// Classical complex symmetric/Hermitian Jacobi eigenbasis; small matrices only.
pub fn real_symmetric_eigen(a:&Matrix)->Result<(Vec<f64>,Matrix),String>{
    if !is_hermitian(a,1e-11) || a.data.iter().any(|z|z.im.abs()>1e-10){
        return Err("symmetric spectral path requires real symmetric matrices".into());
    }
    let n=a.rows;let mut d=a.clone();let mut v=Matrix::identity(n);
    for _ in 0..(120*n*n).max(1) {
        let mut p=0;let mut q=0;let mut max=0.0;
        for i in 0..n {for j in i+1..n {
            let x=d.get(i,j).re.abs();if x>max {max=x;p=i;q=j;}
        }}
        if max<1e-13 {return Ok(((0..n).map(|i|d.get(i,i).re).collect(),v));}
        let app=d.get(p,p).re;let aqq=d.get(q,q).re;let apq=d.get(p,q).re;
        let theta=0.5*(2.0*apq).atan2(aqq-app);
        let c=theta.cos();let s=theta.sin();
        // Rotation columns p,q: v_p <- c*v_p - s*v_q.
        for k in 0..n {
            let vp=v.get(k,p);let vq=v.get(k,q);
            v.set(k,p,vp*c-vq*s);v.set(k,q,vp*s+vq*c);
        }
        for k in 0..n {if k!=p && k!=q {
            let akp=d.get(k,p).re;let akq=d.get(k,q).re;
            let np=c*akp-s*akq;let nq=s*akp+c*akq;
            d.set(k,p,C64::from(np));d.set(p,k,C64::from(np));
            d.set(k,q,C64::from(nq));d.set(q,k,C64::from(nq));
        }}
        d.set(p,p,C64::from(c*c*app-2.0*c*s*apq+s*s*aqq));
        d.set(q,q,C64::from(s*s*app+2.0*c*s*apq+c*c*aqq));
        d.set(p,q,C64::ZERO);d.set(q,p,C64::ZERO);
    }
    Err("Jacobi eigen iteration failed to converge".into())
}
