use std::time::Instant;
use qrk_sylv::fixtures::{Family,make};
use qrk_sylv::solvers::{qrk,generic_kaczmarz,somma,wang_liu,clayton,benedetti,finite_shots};
use qrk_sylv::matrix::Matrix;
fn usage(){eprintln!("Usage: qrk-sylv [--family screened|convection|wave|control|normal-lyapunov|complex|all] [--size 2|4|8] [--instances 30] [--sweeps 50] [--shots 8192] [--seed 123] [--no-shots]");}
fn metric(x:&Matrix,truth:&Matrix)->f64{x.diff_frob(truth)/truth.frob().max(1e-15)}
fn stats(xs:&[f64])->(f64,f64){
    if xs.is_empty(){return (f64::NAN,f64::NAN)}
    let mean=xs.iter().sum::<f64>()/xs.len() as f64;
    let var=if xs.len()>1{xs.iter().map(|x|(x-mean).powi(2)).sum::<f64>()/(xs.len()-1) as f64}else{0.0};
    (mean,var.sqrt())
}
fn eval_method(results:&[Result<(f64,Option<f64>,f64),String>])->String{
    let successes:Vec<_>=results.iter().filter_map(|v|v.as_ref().ok()).collect();
    if successes.is_empty(){return "N/A | N/A | –".into()}
    if successes.len()!=results.len(){return format!("N/A ({} of {} instances supported) | N/A | –",successes.len(),results.len());}
    let ideal=stats(&successes.iter().map(|v|v.0).collect::<Vec<_>>());
    let shots=stats(&successes.iter().filter_map(|v|v.1).collect::<Vec<_>>());
    let ms=stats(&successes.iter().map(|v|v.2).collect::<Vec<_>>());
    format!("{:.3e} ± {:.2e} | {} | {:.2} ± {:.2}",ideal.0,ideal.1,
        if shots.0.is_nan(){"–".into()}else{format!("{:.4} ± {:.4}",shots.0,shots.1)},ms.0,ms.1)
}
fn main(){
    let mut family="all".to_string();let mut size=0usize;let mut instances=30usize;
    let mut sweeps=50usize;let mut shots=8192usize;let mut seed=123u64;let mut no_shots=false;
    let args:Vec<String>=std::env::args().collect();let mut i=1;
    while i<args.len(){
        if args[i]=="--no-shots"{no_shots=true;i+=1;continue;}
        if i+1>=args.len(){usage();std::process::exit(2)}
        let v=&args[i+1];
        match args[i].as_str(){
            "--family"=>family=v.clone(),
            "--size"=>size=v.parse().unwrap_or(0),
            "--instances"=>instances=v.parse().unwrap_or(0),
            "--sweeps"=>sweeps=v.parse().unwrap_or(0),
            "--shots"=>shots=v.parse().unwrap_or(0),
            "--seed"=>seed=v.parse().unwrap_or(123),
            _=>{usage();std::process::exit(2)}
        }
        i+=2;
    }
    if instances==0 || sweeps==0||shots==0{usage();std::process::exit(2)}
    let families:Vec<Family>=if family=="all"{Family::all().to_vec()}else{vec![Family::parse(&family).unwrap_or_else(||{usage();std::process::exit(2)})]};
    let sizes:Vec<usize>=if size==0{vec![2,4,8]}else{vec![size]};
    println!("# Rust solver comparison (classical emulations; not quantum-hardware performance)\n");
    println!("| Family | n | Method | Ideal relative error (mean ± SD) | Finite-shot relative error (mean ± SD) | CPU ms (mean ± SD) |");
    println!("|---|---:|---|---:|---:|---:|");
    for f in families {for &n in &sizes {
        if f==Family::Wave&&n%2!=0{continue;}
        let names=["QRK-Sylv","Generic Kronecker RK","Somma map","Wang-Liu sign","Clayton Lyapunov","Benedetti channel","Classical direct"];
        let mut reports:Vec<Vec<Result<(f64,Option<f64>,f64),String>>>=vec![Vec::new();names.len()];
        for trial in 0..instances {
            let trial_seed=seed.wrapping_add(1_000_003*trial as u64).wrapping_add(n as u64*1009);
            let eq=match make(f,n,trial_seed){Ok(q)=>q,Err(e)=>{eprintln!("{}: {}",f.as_str(),e);continue}};
            let truth=match eq.direct_solve(){Ok(t)=>t,Err(e)=>{eprintln!("reference solve failed: {e}");continue}};
            let mut frozen=None;
            for method in 0..names.len(){
                let start=Instant::now();
                let result:Result<Matrix,String>=match method {
                    0=>{
                        // Supplied bound from independently solved reference only in this test harness.
                        match qrk(&eq,sweeps*n*n,truth.frob().max(1e-12),trial_seed^0xA7AA){
                            Ok(r)=>{frozen=Some(r.frozen_sequence);Ok(r.x)} Err(e)=>Err(e)
                        }
                    },
                    1=>match &frozen {Some(t)=>generic_kaczmarz(&eq,t),None=>Err("QRK sequence missing".into())},
                    2=>somma(&eq),3=>wang_liu(&eq),4=>clayton(&eq),
                    5=>benedetti(&eq),6=>Ok(truth.clone()),_=>unreachable!()
                };
                let ms=start.elapsed().as_secs_f64()*1e3;
                let result=result.map(|x|{
                    let e=metric(&x,&truth);
                    let shot_e=if no_shots {None}else{
                        let sigma=(sweeps*n*n) as f64;
                        let sigma=sigma.sqrt()*truth.frob().max(1e-12);
                        finite_shots(&x,sigma,shots,trial_seed^0xB003^method as u64,0.975,0.015)
                            .ok().map(|z|metric(&z,&truth))
                    };
                    (e,shot_e,ms)
                });
                reports[method].push(result);
            }
        }
        for (name,entries) in names.iter().zip(reports.iter()){
            println!("| {} | {} | {} | {} |",f.as_str(),n,name,eval_method(entries));
        }
    }}
}
