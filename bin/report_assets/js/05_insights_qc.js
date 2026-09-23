// ---- Kraken2: how much of each sample is the organism the cohort is about ----
// target_pct is the share of the CLASSIFIED reads that fall in the cohort's target clade (the one
// most samples are dominated by, e.g. the M. tuberculosis complex). Below KRK_OFF the sample is
// mostly something else; between KRK_OFF and KRK_PURE it carries a real share of another taxon.
var KRK_PURE=90, KRK_OFF=50, KRK_UNCL=15;
// A payload without a target (an older report re-rendered) falls back to the taxon most samples are
// dominated by, measured as its share of the classified reads, which is what target_pct is.
function krkTarget(){var K=R.kraken; if(K&&K.target&&K.target.name)return K.target.name;
  var c={},best=null; ((K&&K.samples)||[]).forEach(function(k){if(k.primary)c[k.primary.name]=(c[k.primary.name]||0)+1;});
  Object.keys(c).forEach(function(n){if(best==null||c[n]>c[best])best=n;}); return best||'the target organism';}
function krkPurity(k){if(k.target_pct!=null)return k.target_pct;
  var t=krkTarget(),hit=null; (k.top||[]).forEach(function(x){if(x.name===t)hit=x;});
  return (hit&&k.classified>0)?Math.min(100,100*hit.pct/k.classified):0;}
function krkState(k){var p=krkPurity(k);return p<KRK_OFF?'off':(p<KRK_PURE?'mixed':(k.unclassified>KRK_UNCL?'uncl':'ok'));}
// the taxa the off-target samples are dominated by, most frequent first: [[name, n], ...]
function krkOffTaxa(off){var c={};off.forEach(function(k){var n=k.primary?k.primary.name:'unclassified';c[n]=(c[n]||0)+1;});
  return Object.keys(c).sort(function(a,b){return c[b]-c[a]||(a<b?-1:1);}).map(function(n){return [n,c[n]];});}
function insKraken(){
  var K=R.kraken; if(!(K&&K.samples&&K.samples.length)) return "";
  var n=K.samples.length, tgt=krkTarget(), off=[], mix=[], uncl=[];
  K.samples.forEach(function(k){var st0=krkState(k); if(st0==='off')off.push(k); else if(st0==='mixed')mix.push(k); else if(st0==='uncl')uncl.push(k);});
  var narr;
  if(!off.length&&!mix.length){
    narr='All <b>'+n+'</b> samples are at least '+KRK_PURE+'% <b>'+esc(tgt)+'</b> (share of the classified reads): no contamination detected by Kraken2.';
  }else{
    narr='';
    if(off.length){var tx=krkOffTaxa(off);
      narr+='<b class="tone-bad">'+off.length+'</b> of '+n+' sample'+(off.length>1?'s are':' is')+' mostly <b>not '+esc(tgt)+'</b>: '+
        tx.slice(0,4).map(function(t){return t[1]+' &#215; <i>'+esc(t[0])+'</i>';}).join(', ')+(tx.length>4?' and '+(tx.length-4)+' other taxa':'')+'.';}
    if(mix.length)narr+=(narr?' ':'')+'<b class="tone-warn">'+mix.length+'</b> '+(off.length?'more ':'')+'carry '+(100-KRK_PURE)+'&#8211;'+(100-KRK_OFF)+'% of other taxa.';
    narr+=' <span class="c">The target is the clade most samples are dominated by; percentages are of the classified reads.</span>';
  }
  if(uncl.length)narr+=' <span class="c">'+uncl.length+' sample'+(uncl.length>1?'s have':' has')+' more than '+KRK_UNCL+'% unclassified reads.</span>';
  var chips=off.concat(mix).map(function(k){var st0=krkState(k);
    return {t:esc(k.s)+' &#183; '+krkPurity(k).toFixed(0)+'%',cls:st0==='off'?'bad':'warn',
            title:krkPurity(k).toFixed(1)+'% '+tgt+(k.primary?'; dominant: '+k.primary.name+' '+k.primary.pct.toFixed(1)+'% of reads':'')};});
  return insBox("Contamination", narr, chips);
}
function insDrug(){
  var d=drSummary(); if(!d) return '';
  var narr,chips=[];
  if(d.nCarriers){
    narr='<b class="tone-bad">'+d.nCarriers+'</b> of '+d.nTyped+' samples carry resistance-associated mutations (WHO 1&#8211;2) that their lineage does not share: '+
      d.drugs.map(function(x){return '<b>'+esc(x.drug)+'</b> in '+x.n;}).join(', ')+'.';
    // MDR: rifampicin and isoniazid both hit in one sample, by mutations beyond the lineage's own
    var rif={},inh={}; d.acquired.forEach(function(m){var ds=(m.dr&&m.dr.length)?m.dr:[m.drug];
      Object.keys(m.samples).forEach(function(s){if(ds.indexOf('RIF')>=0)rif[s]=1;if(ds.indexOf('INH')>=0)inh[s]=1;});});
    var mdr=Object.keys(rif).filter(function(s){return inh[s];});
    if(mdr.length)narr+=' <b class="tone-bad">'+mdr.length+' MDR</b> (rifampicin and isoniazid in one sample).';
    chips=d.acquired.map(function(m){return {t:esc(drMutLab(m))+' &#183; '+m.n,cls:'bad',title:m.drug+', WHO grade '+m.gn+', carried by '+m.n+' samples'};});
  }else{
    narr='<b>No resistance-associated mutation (WHO 1&#8211;2)</b> beyond what each lineage carries, across '+d.nTyped+' samples.';
  }
  if(d.wide.length)narr+=' <span class="c">Shared by a whole lineage and not counted: '+d.wide.map(function(m){return drWideTxt(m,d.nLin);}).join('; ')+'.</span>';
  narr+=' <span class="c">A genomic screen against the WHO catalogue, not a drug-susceptibility result.</span>';
  return insBox('Resistance', narr, chips);
}
// Counts come from the flags, never from the mere presence of a lineage breakdown: every typed sample
// carries one (s.linf), and counting those reported 176 "mixed" samples in a cohort with none.
function insLineages(){
  if(!R.lin_present||!R.samples||!R.samples.length)return '';
  var tot=R.samples.length, groups={}, order=[];
  R.samples.forEach(function(s){if(untyped(s.lineage))return;var k=s.lineage;if(!groups.hasOwnProperty(k)){groups[k]=0;order.push(k);}groups[k]++;});
  order.sort(function(a,b){return groups[b]-groups[a];});
  var typed=order.reduce(function(a,k){return a+groups[k];},0), notTyped=tot-typed;
  var mixed=R.samples.filter(function(s){return s.f.indexOf('MIXED')>=0;});
  var mism=R.samples.filter(function(s){return s.f.indexOf('LINEAGE_MISMATCH')>=0;});
  var narr;
  if(!order.length){
    narr='No sample could be typed to a lineage.';
  }else{
    narr=(order.length===1?'Every typed sample is <b>'+esc(linLabel(order[0]))+'</b>':'Samples type as '+
      order.slice(0,4).map(function(k){return '<b>'+esc(linLabel(k))+'</b> ('+groups[k]+')';}).join(', ')+
      (order.length>4?' and '+(order.length-4)+' more':''))+'.';
    if(notTyped)narr+=' <b class="tone-warn">'+notTyped+'</b> could not be typed.';
    narr+=mixed.length?(' <b class="tone-bad">'+mixed.length+'</b> carr'+(mixed.length>1?'y':'ies')+' a second lineage above '+thr.mixed_min_frac+'% of the markers: a mixed infection or cross-contamination.')
                      :' No sample carries a second lineage.';
    if(mism.length)narr+=' <b class="tone-bad">'+mism.length+'</b> type'+(mism.length>1?'':'s')+' as a different lineage from the other samples on their reference.';
  }
  var chips=mixed.map(function(s){return {t:esc(s.s),cls:'bad',title:flagWhy(s,'MIXED')};})
    .concat(mism.map(function(s){return {t:esc(s.s)+' &#183; '+esc(linMain(s.lineage)),cls:'bad',title:flagWhy(s,'LINEAGE_MISMATCH')};}));
  return insBox('Lineages',narr,chips);
}
function insDist(){
  var dist=R.dist||[];
  if(!dist.length)return '';
  function pct(sorted,p){var n=sorted.length;if(!n)return NaN;var idx=(n-1)*p,lo=Math.floor(idx),hi=Math.ceil(idx);if(lo==hi)return sorted[lo];return sorted[lo]+(sorted[hi]-sorted[lo])*(idx-lo);}
  var scanned=0,withOut=0,best=null,skewBest=null;
  for(var i=0;i<dist.length;i++){
    var pk=dist[i],mt=MET[pk]||{label:pk,kind:'float'};
    var rows=R.samples.filter(function(s){return s.m[pk]!=null&&isFinite(s.m[pk]);});
    if(rows.length<5)continue;
    scanned++;
    var vals=rows.map(function(s){return s.m[pk];});
    var sorted=vals.slice().sort(function(a,b){return a-b;});
    var q1=pct(sorted,0.25),med=_median(vals),q3=pct(sorted,0.75),iqr=q3-q1;
    var bowley=iqr>0?(((q3-med)-(med-q1))/iqr):0;
    var outs=[];
    if(iqr>0){
      var flo=q1-1.5*iqr,fhi=q3+1.5*iqr;
      rows.forEach(function(s){var v=s.m[pk];if(v<flo||v>fhi)outs.push({s:s.s,v:v,vr:s.v,hi:v>fhi});});
    }
    var rec={pk:pk,label:mt.label,kind:mt.kind,med:med,nout:outs.length,outs:outs,bowley:bowley};
    if(outs.length>0)withOut++;
    if(!best||rec.nout>best.nout||(rec.nout==best.nout&&Math.abs(rec.bowley)>Math.abs(best.bowley)))best=rec;
    if(!skewBest||Math.abs(rec.bowley)>Math.abs(skewBest.bowley))skewBest=rec;
  }
  if(!scanned)return '';
  var narr,chips=[];
  if(best&&best.nout>0){
    var os=best.outs.slice().sort(function(a,b){return Math.abs(b.v-best.med)-Math.abs(a.v-best.med);});
    var tone=best.nout>=3?'tone-bad':'tone-warn';
    narr='Across <b>'+scanned+'</b> distributions, <b>'+withOut+'</b> carry at least one sample beyond the 1.5&times;IQR fence. '+
      '<b>'+esc(best.label)+'</b> is the most dispersed: <b class="'+tone+'">'+best.nout+'</b> outlier'+(best.nout>1?'s':'')+
      ' (median '+shortv(best.med,best.kind)+'). Fences flag spread, not QC pass/fail.';
    for(var k=0;k<os.length&&k<6;k++){
      var o=os[k];
      chips.push({t:o.s+' '+shortv(o.v,best.kind),cls:(o.vr=='FAIL'?'bad':'warn'),title:esc(best.label)+' '+shortv(o.v,best.kind)+(o.hi?' (high tail)':' (low tail)')});
    }
    if(os.length>6)chips.push({t:'+'+(os.length-6)+' more',cls:'',title:'additional outlier samples'});
  }else{
    if(!skewBest||Math.abs(skewBest.bowley)<0.05){
      narr='Across <b>'+scanned+'</b> distributions, none breach the 1.5&times;IQR fence and skew is negligible &mdash; the panel is tightly banded.';
    }else{
      var dir=skewBest.bowley>0?'right':'left';
      narr='Across <b>'+scanned+'</b> distributions, none breach the 1.5&times;IQR fence. '+
        '<b>'+esc(skewBest.label)+'</b> shows the strongest skew (<b class="tone-warn">'+(skewBest.bowley>0?'+':'')+skewBest.bowley.toFixed(2)+'</b> Bowley, '+dir+'-tailed; median '+shortv(skewBest.med,skewBest.kind)+').';
    }
  }
  return insBox('Distributions',narr,chips);
}
// Metric pairs that are one quantity measured twice (a percent and its complement, the same coverage
// at two cut-offs, one variant count and its derivatives). Their correlation is arithmetic, so it
// is never the headline: "Mapped % vs Unmapped % at r = -1" told the reader nothing.
var COUPLED=[['mapped_pct','unmapped_pct'],['missing_pct','callable_pct'],['mean_depth','median_depth','est_genome_cov'],
  ['breadth_pct','breadth5x_pct','breadth10x_pct'],['properly_paired_pct','singleton_pct'],['raw_reads','trimmed_reads'],
  ['q20_pct','q30_pct'],
  // every annotation class is a share of the same variant calls, so it grows with the SNP count
  ['snps','total_variants','snp_density','ann_high','ann_moderate','ann_low','ann_modifier']];
function coupled(a,b){for(var i=0;i<COUPLED.length;i++){if(COUPLED[i].indexOf(a)>=0&&COUPLED[i].indexOf(b)>=0)return true;}return false;}
function insScatter(){
  if(!R||!R.samples||R.samples.length<3)return '';
  var keys=[]; (R.dist||[]).forEach(function(k){if(MET[k])keys.push(k);});
  if(keys.length<2)return '';
  function pear(xs,ys){var n=xs.length,i,sx=0,sy=0;for(i=0;i<n;i++){sx+=xs[i];sy+=ys[i];}
    var mx=sx/n,my=sy/n,sxy=0,sxx=0,syy=0;for(i=0;i<n;i++){var dx=xs[i]-mx,dy=ys[i]-my;sxy+=dx*dy;sxx+=dx*dx;syy+=dy*dy;}
    if(sxx<=0||syy<=0)return null;return sxy/Math.sqrt(sxx*syy);}
  var best=null,t,m,xv,yv;
  for(var i=0;i<keys.length;i++){for(var j=i+1;j<keys.length;j++){
    var ka=keys[i],kb=keys[j],xs=[],ys=[];
    if(coupled(ka,kb))continue;
    for(t=0;t<R.samples.length;t++){m=R.samples[t].m;xv=m[ka];yv=m[kb];
      if(xv!=null&&yv!=null&&isFinite(xv)&&isFinite(yv)){xs.push(xv);ys.push(yv);}}
    if(xs.length<5)continue;
    var r=pear(xs,ys);if(r==null)continue;
    if(!best||Math.abs(r)>Math.abs(best.r))best={a:ka,b:kb,r:r,n:xs.length};
  }}
  if(!best)return '';
  var GRP={raw_reads:'seq',trimmed_reads:'seq',
    mean_depth:'cov',median_depth:'cov',breadth_pct:'cov',breadth5x_pct:'cov',breadth10x_pct:'cov',mapped_pct:'cov',unmapped_pct:'cov',properly_paired_pct:'cov',missing_pct:'cov',callable_pct:'cov',coverage_cv:'cov',est_genome_cov:'cov',
    total_variants:'var',snps:'var',het_variants:'var',homo_indels:'var',iupac_pct:'var'};
  var la=MET[best.a].label,lb=MET[best.b].label,ar=best.r,ab=Math.abs(ar);
  var rtxt=(ar>=0?'+':'')+ar.toFixed(2);
  var dir=ar>=0?'positive':'negative';
  var strength=ab>=0.85?'near-perfect':ab>=0.7?'strong':ab>=0.5?'moderate':(ab>=0.3?'weak':'negligible');
  var expected=GRP[best.a]&&GRP[best.b]&&GRP[best.a]===GRP[best.b];
  var surprising=(!expected)&&ab>=0.6;
  var narr='Across '+best.n+' samples, the tightest linear link between two different quantities is <b>'+esc(la)+'</b> vs <b>'+esc(lb)+'</b> ('+
    (surprising?'<b class="tone-warn">Pearson r='+rtxt+'</b>':'<b>Pearson r='+rtxt+'</b>')+', '+strength+' '+dir+'). ';
  if(expected)narr+='These share a QC family and are mechanically coupled, so the correlation is <b>expected</b>.';
  else if(surprising)narr+='<b class="tone-warn">Not an obvious QC pairing</b> — a link worth a look before treating these axes as independent.';
  else narr+='No strong linear structure stands out between the plotted metrics.';
  narr+=' <span style="opacity:.7">Pearson captures only linear trends'+(best.n<8?' and n='+best.n+' is small — read with care':' (n='+best.n+')')+'.</span>';
  var chips=[{t:la,cls:''},{t:lb,cls:''},{t:'r='+rtxt,cls:(surprising?'warn':(expected?'good':'')),title:strength+' '+dir+' correlation over '+best.n+' samples with both metrics'}];
  return insBox('Correlation',narr,chips);
}
function insCorr(){
  if(!R||!R.samples||!R.dist||typeof MET==='undefined') return "";
  var pool=(typeof visible==='function')?visible():R.samples;
  if(!pool||pool.length<4) return "";
  var keys=R.dist.filter(function(k){return MET[k];}).slice(0,16);
  if(keys.length<2) return "";
  var vals={};
  keys.forEach(function(k){vals[k]=pool.map(function(s){return (s&&s.m)?s.m[k]:null;});});
  function rankArr(a){
    var idx=a.map(function(v,i){return i;});
    idx.sort(function(x,y){return a[x]-a[y];});
    var ranks=new Array(a.length),t=0;
    while(t<idx.length){
      var u=t;
      while(u+1<idx.length&&a[idx[u+1]]===a[idx[t]])u++;
      var avg=(t+u)/2+1;
      for(var q=t;q<=u;q++)ranks[idx[q]]=avg;
      t=u+1;
    }
    return ranks;
  }
  function spear(xs,ys){
    var n=xs.length; if(n<4) return null;
    var rx=rankArr(xs),ry=rankArr(ys),mx=0,my=0,i;
    for(i=0;i<n;i++){mx+=rx[i];my+=ry[i];}
    mx/=n; my/=n;
    var num=0,dx=0,dy=0;
    for(i=0;i<n;i++){var a=rx[i]-mx,b=ry[i]-my;num+=a*b;dx+=a*a;dy+=b*b;}
    if(dx<=0||dy<=0) return null;
    return num/Math.sqrt(dx*dy);
  }
  var strong=0,total=0,best=null,bi=-1,bj=-1,bestN=0,i,j,t;
  for(i=0;i<keys.length;i++){
    for(j=i+1;j<keys.length;j++){
      if(coupled(keys[i],keys[j]))continue;   // one quantity measured twice: its rho is arithmetic, not a finding
      var xs=[],ys=[],vi=vals[keys[i]],vj=vals[keys[j]];
      for(t=0;t<pool.length;t++){var xv=vj[t],yv=vi[t]; if(xv!=null&&yv!=null){xs.push(xv);ys.push(yv);}}
      var r=spear(xs,ys);
      if(r==null) continue;
      total++;
      if(Math.abs(r)>=0.7) strong++;
      if(best==null||Math.abs(r)>Math.abs(best)){best=r;bi=i;bj=j;bestN=xs.length;}
    }
  }
  if(best==null||total===0) return "";
  var labA=esc(MET[keys[bi]].label),labB=esc(MET[keys[bj]].label);
  var absr=Math.abs(best).toFixed(2);
  var rhoStr=(best>=0?'+':'−')+absr;
  var redundant=Math.abs(best)>=0.9;
  var narr;
  if(strong>0){
    narr='Across '+pool.length+' samples, <b'+(redundant?' class="tone-warn"':'')+'>'+strong+' of '+total+'</b> pairs of different metrics are tightly coupled (|ρ| ≥ 0.7; a percent and its complement are left out). Strongest: <b>'+labA+'</b> vs <b>'+labB+'</b> at ρ '+rhoStr+
      (redundant?' — near-redundant.':'.')+' Spearman rank over '+bestN+' shared points, so with few samples high |ρ| can arise by chance.';
  } else {
    narr='None of '+total+' metric pairs reach |ρ| ≥ 0.7 across '+pool.length+' samples — QC metrics vary largely independently. Strongest is <b>'+labA+'</b> vs <b>'+labB+'</b> at ρ '+rhoStr+' (Spearman rank, '+bestN+' shared points).';
  }
  var chips=[];
  chips.push({t:MET[keys[bi]].label+' ~ '+MET[keys[bj]].label+'  ρ '+rhoStr, cls:(Math.abs(best)>=0.7?(redundant?'bad':'warn'):''), title:'Spearman ρ over '+bestN+' shared samples'});
  chips.push({t:strong+' pair'+(strong===1?'':'s')+' |ρ|≥0.7', cls:(strong>0?'warn':'good'), title:strong+' of '+total+' evaluated metric pairs cross the |ρ|≥0.7 threshold'});
  return insBox('Correlation', narr, chips);
}
function insQCspace(){
  if(!R||!R.samples||!R.samples.length)return '';
  var cut=(R.mahal_cut!=null)?R.mahal_cut:null;
  var rows=[],i;
  for(i=0;i<R.samples.length;i++){var s=R.samples[i];if(s&&s.m&&s.m.qc_mahal!=null)rows.push(s);}
  if(!rows.length)return '';
  rows.sort(function(a,b){return b.m.qc_mahal-a.m.qc_mahal;});
  var out=[];
  if(cut!=null){for(i=0;i<rows.length;i++){if(rows[i].m.qc_mahal>cut)out.push(rows[i]);}}
  var top=rows[0];
  var pev1=null;
  try{if(typeof pca2==='function'&&typeof visible==='function'){var P=pca2(visible());if(P&&P.ok&&P.pev&&P.pev.length)pev1=P.pev[0];}}catch(e){}
  var drvTxt='';
  if(top._mdrv&&top._mdrv.length){var d0=top._mdrv[0];drvTxt=', driven by <b>'+esc(d0.label)+'</b> ('+(d0.z>=0?'+':'')+d0.z.toFixed(1)+'σ)';}
  var pevTxt=(pev1!=null)?' PC1 captures <b>'+pev1.toFixed(1)+'%</b> of the QC-metric variance.':'';
  var narr,chips=[];
  if(cut!=null&&out.length){
    narr='<b class="tone-warn">'+out.length+'</b> of '+rows.length+' samples sit <b>outside</b> the QC-space envelope (Mahalanobis d² &gt; '+cut.toFixed(1)+' cut). <b>'+esc(top.s)+'</b> is the most extreme at d² '+shortv(top.m.qc_mahal,'float')+drvTxt+'.'+pevTxt;
    out.forEach(function(s){chips.push({t:s.s+' · d² '+shortv(s.m.qc_mahal,'float'),cls:(s.v=='FAIL'?'bad':'warn'),title:'Mahalanobis distance in QC-metric space exceeds the '+cut.toFixed(1)+' cut'});});
  }else if(cut!=null){
    narr='All <b>'+rows.length+'</b> samples fall <b>within</b> the QC-space envelope (d² ≤ '+cut.toFixed(1)+' cut). <b>'+esc(top.s)+'</b> lies closest to the edge at d² '+shortv(top.m.qc_mahal,'float')+drvTxt+'.'+pevTxt;
    chips.push({t:top.s+' · d² '+shortv(top.m.qc_mahal,'float'),cls:'',title:'nearest to the QC-space boundary'});
  }else{
    narr='<b>'+esc(top.s)+'</b> is the most unusual sample in QC-metric space, at Mahalanobis d² '+shortv(top.m.qc_mahal,'float')+drvTxt+' (no cut defined — ranking only).'+pevTxt;
    chips.push({t:top.s+' · d² '+shortv(top.m.qc_mahal,'float'),cls:'warn',title:'highest Mahalanobis distance in view'});
  }
  return insBox('Outliers',narr,chips);
}
function insRefBias(){
  if(typeof R==='undefined'||!R||!R.samples||!R.samples.length)return '';
  var thr=R.thresholds||{};var gx=thr.missing_max;
  if(gx==null)return '';
  var yk=R.snp_density_ok?'snp_density':'snps';
  var ykind=R.snp_density_ok?'float':'int';
  var ylab=R.snp_density_ok?'SNPs/callable Mb':'SNPs vs ref';
  function xof(s){var m=s.m||{};return m.missing_pct!=null?m.missing_pct:(m.callable_pct!=null?100-m.callable_pct:null);}
  function pct(v){return (+v).toFixed(1)+'%';}
  var pool=[],i,s;
  for(i=0;i<R.samples.length;i++){s=R.samples[i];if(s&&s.m&&xof(s)!=null&&s.m[yk]!=null)pool.push(s);}
  if(pool.length<3)return '';
  var ys=[];for(i=0;i<pool.length;i++)ys.push(pool[i].m[yk]);
  var scr=divScreen(ys),dmed=scr.med,ylo=scr.lo;
  if(!scr.on){
    var far=pool.filter(function(x){return x.f.indexOf('SNP_HIGH')>=0;});
    return insBox('Divergence','The cohort sits almost on its reference (median <b>'+fmt(dmed,ykind)+' '+ylab+'</b>), as when samples are mapped to their own ancestor. '+
      'Against such a reference a sample with fewer SNPs than the rest is noise, so no reference-bias screen is run.'+
      (far.length?' <b class="tone-warn">'+far.length+'</b> sample'+(far.length>1?'s sit':' sits')+' far above the rest instead (<i>Unusually many SNPs</i>): check that they were mapped to the right genome.':''),
      far.map(function(x){return {t:esc(x.s)+' &#183; '+fmt(x.m[yk],ykind),cls:'warn',title:fmt(x.m[yk],ykind)+' '+ylab};}));
  }
  var refbias=[],lowcov=[],xv,yv;
  for(i=0;i<pool.length;i++){s=pool[i];xv=xof(s);yv=s.m[yk];
    if(yv<ylo){if(xv<=gx)refbias.push(s);else lowcov.push(s);}}
  var medTxt=fmt(dmed,ykind)+' '+ylab;var chips,narr;
  if(refbias.length){
    chips=refbias.map(function(s){return {t:esc(s.s),cls:'bad',title:fmt(s.m[yk],ykind)+' '+ylab+' at '+pct(xof(s))+' missing'};});
    narr='<b class="tone-bad">'+refbias.length+'</b> sample'+(refbias.length>1?'s':'')+' sit'+(refbias.length>1?'':'s')+' anomalously low on divergence (&lt;'+fmt(ylo,ykind)+' vs '+medTxt+' cohort median) yet map near-completely (&le;'+pct(gx)+' missing) &mdash; the complete-but-quiet signature of <b>reference bias</b> or reference-ward contamination rather than honest low coverage. Robust MAD cut; small cohorts are noisy, so confirm on the plot.';
    return insBox('Divergence',narr,chips);
  }
  if(lowcov.length){
    chips=lowcov.map(function(s){return {t:esc(s.s),cls:'warn',title:fmt(s.m[yk],ykind)+' '+ylab+' at '+pct(xof(s))+' missing'};});
    narr='No reference-bias suspects: the '+(lowcov.length>1?lowcov.length+' low-divergence outliers all carry':'lone low-divergence outlier carries')+' elevated missingness (&gt;'+pct(gx)+'), so the sparse calls are honestly explained by <b class="tone-warn">incomplete coverage</b>, not bias (median '+medTxt+').';
    return insBox('Divergence',narr,chips);
  }
  narr='All <b>'+pool.length+'</b> samples track the expected divergence&ndash;completeness trend (median '+medTxt+'); none show the complete-but-low-divergence pattern that flags <b>reference bias</b> or contamination.';
  return insBox('Divergence',narr);
}
