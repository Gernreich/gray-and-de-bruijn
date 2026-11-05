// Count every least-turning grid ring: python3 walks.py --graph | (rustc -O -o /tmp/cnt tools/count_grid_rings.rs && /tmp/cnt)
// Prints directed and undirected counts (2161152 undirected; about 13 minutes).
use std::io::Read;
struct G { n: usize, adj: Vec<Vec<(usize,u32)>>, start: usize, bound: u32 }
fn ok(g:&G, seen:u32, cur:usize) -> bool {
  // every unvisited vertex needs 2 usable neighbours (unvisited, or cur, or start)
  for v in 0..g.n { if seen>>v&1==1 {continue}
    let mut d=0; for &(w,_) in &g.adj[v] { if seen>>w&1==0 || w==cur || w==g.start { d+=1; } }
    if d<2 {return false} }
  true
}
fn dfs(g:&G, v:usize, seen:u32, n:usize, cost:u32, cnt:&mut u64) {
  if cost + (g.n-n) as u32 > g.bound {return}
  if n==g.n { for &(w,c) in &g.adj[v] { if w==g.start && cost+c==g.bound { *cnt+=1; } } return }
  if !ok(g,seen,v) {return}
  for &(w,c) in &g.adj[v] { if seen>>w&1==0 { dfs(g,w,seen|1<<w,n+1,cost+c,cnt); } }
}
fn main(){
  let mut s=String::new(); std::io::stdin().read_to_string(&mut s).unwrap();
  let mut l=s.lines(); let h:Vec<usize>=l.next().unwrap().split(' ').map(|x|x.parse().unwrap()).collect();
  let adj=l.map(|x| x.split(' ').filter(|t|!t.is_empty()).map(|t|{let mut p=t.split(':');(p.next().unwrap().parse().unwrap(),p.next().unwrap().parse().unwrap())}).collect()).collect();
  let g=G{n:h[0],adj,start:h[1],bound:32};
  let mut cnt=0; dfs(&g,g.start,1<<g.start,1,0,&mut cnt);
  println!("directed {} undirected {}",cnt,cnt/2);
}
