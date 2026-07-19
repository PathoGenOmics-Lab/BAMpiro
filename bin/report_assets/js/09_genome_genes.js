var GTRACKS=[{k:'missing',lab:'Missing',base:[214,64,58]},{k:'snp',lab:'SNPs',base:[14,139,168]},
             {k:'het',lab:'Het',base:[221,138,26]},{k:'indel',lab:'Indels',base:[124,92,191]}];
function gtrackHas(k){return R.samples.some(function(s){return k=='missing'?!!s.miss:(s.trk&&s.trk[k]);});}
function gtrackMax(k){var mx=0;R.samples.forEach(function(s){var p=k=='missing'?s.miss:(s.trk&&s.trk[k]);if(p)p.forEach(function(v){if(v!=null&&v>mx)mx=v;});});return mx;}
function gcol(base,mv){if(mv==null)return TH.cellnull;if(mv>1)mv=1;var a=isDark()?[30,42,56]:[238,244,240];
  return 'rgb('+Math.round(a[0]+(base[0]-a[0])*mv)+','+Math.round(a[1]+(base[1]-a[1])*mv)+','+Math.round(a[2]+(base[2]-a[2])*mv)+')';}
// ---- Functional annotation (snpEff impact + effect classes per sample) ----
var IMP=[['ann_high','HIGH','#e0544f'],['ann_moderate','MODERATE','#e6b25a'],
         ['ann_low','LOW','#4f83c2'],['ann_modifier','MODIFIER','#cbd5e1']];
var EFFCLS=[['eff_missense','missense','#4f83c2'],['eff_synonymous','synonymous','#22a06b'],
            ['eff_stop_gained','stop gained','#e0544f'],['eff_frameshift','frameshift','#c0392b'],
            ['eff_start_lost','start lost','#b5651d'],['eff_stop_lost','stop lost','#d98c2b'],
            ['eff_inframe_indel','inframe indel','#7c5cbf'],['eff_splice','splice','#a0508a'],
            ['eff_intergenic','intergenic','#9aa7b6'],['eff_regulatory','regulatory','#c4ccd7']];
var FUNCTION_CAPTION='snpEff functional annotation per sample. The bar splits each sample&#39;s annotated variants into snpEff impact classes (HIGH, MODERATE, LOW, MODIFIER; most severe first). Impact is snpEff&#39;s own severity call from the reference gene model, so these counts inherit any error in the reference annotation, and a HIGH-impact excess can be genuine loss-of-function, a wrong reference database, contamination, or an indel-calling artifact rather than biology. The missense/silent ratio shown is a pN/pS PROXY (a spectrum sanity check within a lineage), NOT dN/dS and NOT a test of selection; it is entangled with divergence and with how snpEff bins effects, so read it alongside the Ti/Tv ratio and use the eskaks pN/pS panel for a real estimate. Counts are absolute (a divergent lineage has more of every class); the excess-HIGH and excess-LoF flags are cohort-relative robust-z, not absolute cut-offs.';
function renderFunction(){
  var host=el('fn_stacks'),cap=el('fn_caption'),leg=el('fn_legend'),coh=el('fn_cohort');
  if(!host)return;
  var rows=R.samples.filter(function(s){return IMP.some(function(p){return s.m[p[0]]!=null;});})
    .slice().sort(function(a,b){return (b.m.ann_high||0)-(a.m.ann_high||0);});
  if(!rows.length){host.innerHTML='<div class="pad nd">no snpEff annotation in the summary (pass annotated VCFs upstream).</div>';if(cap)cap.innerHTML='';if(leg)leg.innerHTML='';if(coh)coh.innerHTML='';return;}
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  host.innerHTML=rows.map(function(s){
    var tot=IMP.reduce(function(a,p){return a+(s.m[p[0]]||0);},0)||1;
    var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s];
    var segs=IMP.map(function(p){var v=s.m[p[0]]||0,w=100*v/tot;if(w<=0)return '';
      return '<div style="width:'+w.toFixed(2)+'%;background:'+p[2]+'" title="'+p[1]+': '+Math.round(v).toLocaleString('en-US')+'"></div>';}).join('');
    var hi=s.m.ann_high!=null?Math.round(s.m.ann_high).toLocaleString('en-US'):'-';
    return '<div class="stack" data-s="'+esc(s.s)+'" style="opacity:'+(dim?0.25:1)+(st.hi==s.s?';background:'+TH.hl:'')+'">'+
      '<span class="sl">'+esc(s.s)+'</span><div class="sb">'+segs+'</div>'+
      '<span class="sv" title="HIGH-impact count">'+hi+'</span></div>';}).join('');
  Array.prototype.forEach.call(host.querySelectorAll('.stack'),function(d){d.onclick=function(){setHi(d.getAttribute('data-s'));};});
  if(leg)leg.innerHTML=IMP.map(function(p){return '<span><i style="background:'+p[2]+'"></i>'+p[1]+'</span>';}).join('')+'<span style="margin-left:auto">right number = HIGH-impact count; sorted worst first</span>';
  var pool=visible().filter(function(s){return EFFCLS.some(function(p){return s.m[p[0]]!=null;});});
  var sums={},grand=0;EFFCLS.forEach(function(p){var t=0;pool.forEach(function(s){if(s.m[p[0]]!=null)t+=s.m[p[0]];});sums[p[0]]=t;grand+=t;});
  var mx=1;EFFCLS.forEach(function(p){if(sums[p[0]]>mx)mx=sums[p[0]];});
  var effRows=EFFCLS.filter(function(p){return sums[p[0]]>0;}).map(function(p){var v=sums[p[0]];
    return '<div class="drow"><div class="dk">'+esc(p[1])+'</div><div class="dbarwrap"><div class="dbar" style="width:'+(100*v/mx).toFixed(1)+'%;background:'+p[2]+'"></div></div>'+
      '<div class="dv">'+Math.round(v).toLocaleString('en-US')+'</div><div class="dp">'+(grand?Math.round(100*v/grand):0)+'%</div></div>';}).join('');
  var mis=sums.eff_missense||0,syn=sums.eff_synonymous||0,ratio=(syn>0)?(mis/syn):null;
  var strip='<div class="cur-read" style="margin:6px 0 10px"><b>'+(ratio==null?'NA':ratio.toFixed(2))+'</b> cohort missense/silent ratio '+
    '<span style="color:#94a3b8;font-weight:400">(pN/pS proxy; '+Math.round(mis).toLocaleString('en-US')+' missense / '+Math.round(syn).toLocaleString('en-US')+' synonymous; not a selection test)</span></div>';
  if(coh)coh.innerHTML='<div style="flex-basis:100%"><div class="dsub" style="margin:2px 0 4px">Cohort effect classes <span style="font-weight:400;color:#94a3b8">(summed over samples in view)</span></div>'+strip+effRows+'</div>';
  if(cap)cap.innerHTML=FUNCTION_CAPTION;
  host.style.setProperty('--stackh',fnZoom+'px'); host.style.setProperty('--sbh',Math.round(fnZoom*0.5)+'px');
  var fz=el('fnzoom'); if(fz){ fz.value=fnZoom; fz.oninput=function(){ fnZoom=+this.value; host.style.setProperty('--stackh',fnZoom+'px'); host.style.setProperty('--sbh',Math.round(fnZoom*0.5)+'px'); }; }
}
var fnZoom=26, genomeZoom=0;   // function per-sample bar row height (px, CSS var); genome row height override (0 = auto)
function renderGenome(){
  var host=el('genome_body'); if(!host)return;
  var samp=R.samples.filter(function(s){return s.miss||s.trk;});
  if(!samp.length){host.innerHTML='<span class="nd" style="padding:0">no consensus/variant data for a genome landscape.</span>';return;}
  var tk=st.gtrack||'missing'; if(!gtrackHas(tk)){ var av=GTRACKS.filter(function(g){return gtrackHas(g.k);}); tk=av.length?av[0].k:'missing'; st.gtrack=tk; }   // fall back to the first track that HAS data (e.g. SNP density when there's no consensus/missing track) instead of showing a spurious 'no data'
  var meta=GTRACKS[0]; GTRACKS.forEach(function(g){if(g.k==tk)meta=g;}); var base=meta.base;
  var rows0=samp.filter(function(s){return (tk=='missing')?s.miss:(s.trk&&s.trk[tk]);});
  if(!rows0.length){host.innerHTML='<span class="nd" style="padding:0">no data for this track.</span>';return;}
  var nb=(tk=='missing')?rows0[0].miss.length:rows0[0].trk[tk].length, gl=R.genome_len||nb;
  var mx=(tk=='missing')?1:(gtrackMax(tk)||1);
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});
  var rows=rows0.slice().sort(function(a,b){var la=(R.lineages||[]).indexOf(a.lineage),lb=(R.lineages||[]).indexOf(b.lineage);
    if(la<0)la=999; if(lb<0)lb=999; if(la!=lb)return la-lb; return a.s.localeCompare(b.s);});
  var W=Math.max(280,host.clientWidth||900), gut=54, profH=54, rowH=genomeZoom>0?genomeZoom:Math.max(3,Math.min(9,Math.floor(300/rows.length))), plotW=W-gut-8;
  var z0=0,z1=nb-1; if(st.gzoom){z0=Math.max(0,Math.min(nb-1,st.gzoom.b0|0));z1=Math.max(z0,Math.min(nb-1,st.gzoom.b1|0));}
  var winN=z1-z0+1;                                            // zoom window (bins); the whole reference when st.gzoom is null
  function x(i){return gut+(i-z0)/winN*plotW;} var bw=plotW/winN+0.6;
  function raw(s,i){if(tk=='missing')return s.miss?s.miss[i]:null;var p=s.trk&&s.trk[tk];return p?p[i]:null;}
  var visRows=rows.filter(function(s){return vis[s.s];}); if(!visRows.length)visRows=rows;
  var agg=[];for(var i=0;i<nb;i++){var sum=0,c=0;visRows.forEach(function(s){var vv=raw(s,i);if(vv!=null){sum+=vv;c++;}});agg.push(c?sum/c:null);}
  var mb=(st.maskOn&&R.mask_bins)?R.mask_bins:null;
  if(mb)for(var mi=0;mi<nb;mi++)if(mb[mi]>=0.5)agg[mi]=null;   // masked bins drop out of the density profile
  var y0=profH+14, totH=y0+rows.length*rowH+18;
  // cohort DENSITY PROFILE (area + line) of the selected track along the reference
  var pmax=Math.max.apply(null,agg.slice(z0,z1+1).map(function(v){return v||0;}).concat([1e-9]));
  var pts=[];for(var pi=z0;pi<=z1;pi++){var pv=agg[pi]||0;pts.push(x(pi).toFixed(1)+','+(profH-(pv/pmax)*(profH-10)).toFixed(1));}
  var bs='rgb('+base[0]+','+base[1]+','+base[2]+')',bf='rgba('+base[0]+','+base[1]+','+base[2]+',.15)';
  var TLAB={missing:'missing %',snp:'SNP density',het:'het density',indel:'indel density'};
  var svg='<svg width="'+W+'" height="'+totH+'">';
  svg+='<line x1="'+gut+'" y1="'+profH+'" x2="'+(gut+plotW).toFixed(1)+'" y2="'+profH+'" stroke="'+TH.grid+'"/>';
  svg+='<path d="M '+gut.toFixed(1)+' '+profH+' L '+pts.join(' L ')+' L '+(gut+plotW).toFixed(1)+' '+profH+' Z" fill="'+bf+'"/>';
  svg+='<polyline points="'+pts.join(' ')+'" fill="none" stroke="'+bs+'" stroke-width="1.3" stroke-linejoin="round"/>';
  svg+='<text x="0" y="11" font-size="9" font-weight="600" fill="'+TH.mut+'">'+TLAB[tk]+'</text>';
  svg+='<text x="'+(gut+plotW).toFixed(1)+'" y="11" font-size="8.5" fill="#94a3b8" text-anchor="end">peak '+(tk=='missing'?(pmax*100).toFixed(0)+'%':(Math.round(pmax*10)/10)+'/bin')+'</text>';
  for(var pj=z0;pj<=z1;pj++)svg+='<rect x="'+x(pj).toFixed(1)+'" y="0" width="'+bw.toFixed(1)+'" height="'+profH+'" fill="transparent" data-bin="'+pj+'" data-v="'+(agg[pj]==null?'':(tk=='missing'?(agg[pj]*100).toFixed(0):(Math.round(agg[pj]*10)/10)))+'"/>';
  rows.forEach(function(s,r){var dim=(st.q||st.onlyFlagged||st.flagFilter||st.ancOnly||st.linFilter)&&!vis[s.s],yy=y0+r*rowH;
    svg+='<rect x="'+(gut-9)+'" y="'+yy+'" width="5" height="'+(rowH-0.4).toFixed(1)+'" fill="'+linColor(s.lineage)+'"'+(dim?' opacity="0.25"':'')+'/>';
    for(var i=z0;i<=z1;i++){var vv=raw(s,i),mv=vv==null?null:(tk=='missing'?vv:vv/mx);
      svg+='<rect x="'+x(i).toFixed(1)+'" y="'+yy+'" width="'+bw.toFixed(1)+'" height="'+(rowH-0.4).toFixed(1)+'" fill="'+gcol(base,mv)+'"'+(dim?' opacity="0.3"':'')+' data-s="'+esc(s.s)+'" data-bin="'+i+'" data-v="'+(vv==null?'':(tk=='missing'?(vv*100).toFixed(0):vv))+'"/>';}
  });
  var yb=y0+rows.length*rowH;
  if(mb)for(var mk=z0;mk<=z1;mk++){var mf=mb[mk]||0;if(mf<=0)continue;   // grey out the masked zones
    svg+='<rect x="'+x(mk).toFixed(1)+'" y="0" width="'+bw.toFixed(1)+'" height="'+yb+'" fill="'+TH.faint+'" opacity="'+(0.12+0.42*mf).toFixed(2)+'"/>';}
  [0,0.25,0.5,0.75,1].forEach(function(t){var px=gut+t*plotW,pos=Math.round((z0+t*winN)/nb*gl);
    svg+='<line x1="'+px.toFixed(1)+'" y1="'+yb+'" x2="'+px.toFixed(1)+'" y2="'+(yb+4)+'" stroke="'+TH.axis+'"/>'+
      '<text x="'+px.toFixed(1)+'" y="'+(yb+14)+'" font-size="8.5" fill="#94a3b8" text-anchor="'+(t==0?'start':t==1?'end':'middle')+'">'+fmtpos(pos)+'</text>';});
  if(tk=='snp'&&st.geneMark&&st.geneMark.b1>=z0&&st.geneMark.b0<=z1){var gm=st.geneMark,mx0=Math.max(gut,x(gm.b0)),mx1=Math.min(gut+plotW,x(gm.b1+1));
    svg+='<rect x="'+mx0.toFixed(1)+'" y="0" width="'+Math.max(2,mx1-mx0).toFixed(1)+'" height="'+yb+'" fill="none" stroke="'+TH.ink+'" stroke-width="1.2" stroke-dasharray="3 2"/>'+
      '<text x="'+Math.min(W-2,(mx0+mx1)/2).toFixed(1)+'" y="'+(profH+11)+'" font-size="9" font-weight="600" fill="'+TH.ink+'" text-anchor="middle">'+esc(gm.name||'')+'</text>';}
  svg+='<rect id="gbrush" x="0" y="0" width="0" height="'+yb+'" fill="rgba(14,139,168,.12)" stroke="#0e8ba8" stroke-width="1" stroke-dasharray="3 2" pointer-events="none" style="display:none"/>';
  GGEO={gut:gut,plotW:plotW,nb:nb,gl:gl,yb:yb,z0:z0,winN:winN};
  host.innerHTML=svg+'</svg>';
  // selection overlay only when it is a SUB-region of the current view (when zoomed to exactly the selection, the whole plot is it)
  if(st.gsel&&!(st.gsel.b0<=z0&&st.gsel.b1>=z1)){var gbx=el('gbrush');if(gbx){var bx0=Math.max(gut,x(st.gsel.b0)),bx1=Math.min(gut+plotW,x(st.gsel.b1+1));
    if(bx1>bx0){gbx.setAttribute('x',bx0.toFixed(1));gbx.setAttribute('width',(bx1-bx0).toFixed(1));gbx.style.display='';}else gbx.style.display='none';}}
  var gz=el('genzoom'); if(gz){ if(genomeZoom<=0)gz.value=rowH; gz.oninput=function(){ genomeZoom=+this.value; renderGenome(); }; }   // row-height zoom (re-renders; slider lives in the h2 so the drag survives)
}

// ---- SNP-dense gene / region detection (cohort SNP density along the reference) ----
function cohortSnp(){var nb=R.nbins||200,agg=zeros(nb),has=false,mb=(st.maskOn&&R.mask_bins)?R.mask_bins:null;
  R.samples.forEach(function(s){var p=s.trk&&s.trk.snp;if(p){has=true;for(var i=0;i<nb&&i<p.length;i++){if(mb&&mb[i]>=0.5)continue;agg[i]+=p[i]||0;}}});
  return has?agg:null;}
function robustMS(a){var v=a.slice().sort(function(x,y){return x-y;}),n=v.length;if(!n)return[0,1];
  var med=n%2?v[(n-1)/2]:(v[n/2-1]+v[n/2])/2,d=a.map(function(x){return Math.abs(x-med);}).sort(function(x,y){return x-y;});
  var mad=n%2?d[(n-1)/2]:(d[n/2-1]+d[n/2])/2,sig=1.4826*mad;
  if(sig<=0){var m=a.reduce(function(p,q){return p+q;},0)/n;sig=Math.sqrt(a.reduce(function(p,q){return p+(q-m)*(q-m);},0)/n)||1;}
  return[med,sig];}
function renderHotspots(){
  var host=el('hot_body'),tbl=el('hottable'),note=el('hot_note'); if(!host||!tbl)return;
  var agg=cohortSnp();
  if(!agg){host.style.display='none';return;} host.style.display='';
  var nb=R.nbins||agg.length, gl=R.genome_len||nb, binbp=gl/nb, genes=R.genes||[], head, body;
  function b0of(p){return Math.max(0,Math.min(nb-1,Math.floor(p/gl*nb)));}
  function geneMaskedFrac(s,e){var iv=R.mask_iv;if(!iv)return 0;var len=(e-s+1)||1,ov=0;
    for(var i=0;i<iv.length;i++){var a=iv[i][0],b=iv[i][1];if(b<=s||a>=e)continue;ov+=Math.max(0,Math.min(e,b)-Math.max(s,a));}
    return ov/len;}
  if(genes.length){
    var items=genes.map(function(g){var s=Math.max(1,g.start),e=Math.max(s,g.end),len=(e-s+1)||1,b0=b0of(s),b1=b0of(e),tot=0;
      for(var b=b0;b<=b1;b++){var binS=b*binbp,binE=(b+1)*binbp,ov=Math.max(0,Math.min(e,binE)-Math.max(s,binS));tot+=agg[b]*(binbp>0?ov/binbp:1);}
      return {name:g.name,start:s,end:e,len:len,snp:tot,dens:tot/(len/1000),b0:b0,b1:b1};});
    if(st.maskOn&&R.mask_iv)items=items.filter(function(x){return geneMaskedFrac(x.start,x.end)<0.5;});   // drop masked genes
    if(st.gsel)items=items.filter(function(x){return x.b1>=st.gsel.b0&&x.b0<=st.gsel.b1;});   // brushed span only
    if(st.hotq){var _hq=st.hotq.toLowerCase();items=items.filter(function(x){return (x.name||'').toLowerCase().indexOf(_hq)>=0;});}   // gene search
    var ms=robustMS(items.map(function(x){return x.dens;}));
    items.forEach(function(x){x.z=ms[1]?(x.dens-ms[0])/ms[1]:0;});
    items.sort(function(a,b){return b.dens-a.dens;});
    var _hcap=st.hotq?300:((el('hotspotsPanel')&&el('hotspotsPanel').classList.contains('expanded'))?150:15);
    head='<tr><th class="s" style="text-align:left">Gene</th><th style="text-align:left">Position</th><th>Length</th><th>Cohort SNPs</th><th>SNPs/kb</th><th>robust z</th></tr>';
    body=items.slice(0,_hcap).map(function(x){return '<tr data-b0="'+x.b0+'" data-b1="'+x.b1+'" data-name="'+esc(x.name)+'">'+
      '<td class="s" style="text-align:left">'+esc(x.name)+geneRvTag(x.name)+'</td><td style="text-align:left">'+fmtpos(x.start)+' - '+fmtpos(x.end)+'</td>'+
      '<td>'+Math.round(x.len).toLocaleString('en-US')+'</td><td>'+Math.round(x.snp).toLocaleString('en-US')+'</td><td>'+x.dens.toFixed(1)+'</td>'+
      '<td'+(x.z>=3?' style="color:var(--fail);font-weight:600"':'')+'>'+(x.z>0?'+':'')+x.z.toFixed(1)+'</td></tr>';}).join('');
    note.textContent=(st.hotq?('Showing '+Math.min(items.length,_hcap)+' of '+items.length+' matching genes. '):'')+'Genes ranked by cohort SNP density (SNPs per kb, summed over all samples); z = robust outlier vs the gene distribution. Resolution is limited to the '+(binbp/1000).toFixed(0)+' kb profile bins.';
    if(!items.length)body='<tr><td class="na" colspan="6" style="text-align:left;padding:8px">no gene matches "'+esc(st.hotq||'')+'".</td></tr>';
  }else{
    var ms2=robustMS(agg),thr=ms2[0]+3*ms2[1],hot=[];
    for(var i=0;i<nb;i++)if(agg[i]>thr){if(hot.length&&hot[hot.length-1].b1==i-1)hot[hot.length-1].b1=i;else hot.push({b0:i,b1:i});}
    hot.forEach(function(h){h.snp=0;for(var b=h.b0;b<=h.b1;b++)h.snp+=agg[b];});
    if(st.gsel)hot=hot.filter(function(h){return h.b1>=st.gsel.b0&&h.b0<=st.gsel.b1;});
    hot.sort(function(a,b){return b.snp-a.snp;});
    head='<tr><th class="s" style="text-align:left">Region</th><th>Bins</th><th>Cohort SNPs</th></tr>';
    body=hot.slice(0,15).map(function(h){return '<tr data-b0="'+h.b0+'" data-b1="'+h.b1+'"><td class="s" style="text-align:left">'+fmtpos(Math.round(h.b0*binbp))+' - '+fmtpos(Math.round((h.b1+1)*binbp))+'</td><td>'+(h.b1-h.b0+1)+'</td><td>'+Math.round(h.snp).toLocaleString('en-US')+'</td></tr>';}).join('');
    note.textContent='SNP-dense regions (robust-z > 3 over the '+(binbp/1000).toFixed(0)+' kb bins). Pass a gene GFF (--gff) to name the genes.';
  }
  if(st.maskOn&&R.mask_bins)note.textContent+=' Masked regions excluded.';
  if(st.gsel)note.textContent+=' Restricted to the brushed span.';
  tbl.innerHTML='<thead>'+head+'</thead><tbody>'+(body||'<tr><td colspan="'+(genes.length?6:3)+'" style="text-align:left;color:#16a34a;padding:8px">no SNP-dense outliers.</td></tr>')+'</tbody>';
  Array.prototype.forEach.call(tbl.querySelectorAll('tbody tr[data-b0]'),function(tr){tr.onclick=function(){
    var _gb0=+tr.getAttribute('data-b0'),_gb1=+tr.getAttribute('data-b1'),_nb=R.nbins||200,_gp=Math.max(3,Math.round((_gb1-_gb0)*0.6)+2);
    st.gtrack='snp'; st.geneMark={b0:_gb0,b1:_gb1,name:tr.getAttribute('data-name')||''}; st.gzoom={b0:Math.max(0,_gb0-_gp),b1:Math.min(_nb-1,_gb1+_gp)};   // zoom to the clicked gene
    Array.prototype.forEach.call(document.querySelectorAll('#gtrack button'),function(x){x.classList.toggle('on',x.getAttribute('data-gt')=='snp');});
    renderGenome(); el('genome').scrollIntoView();};});
}

// ---- per-lineage summary + cohort composition bar ----
function med(a){if(!a.length)return null;var s=a.slice().sort(function(x,y){return x-y;});return s[Math.floor(s.length/2)];}
function renderLineages(){
  if(!R.lin_present)return;
  var groups={}; R.samples.forEach(function(s){var k=s.lineage||'NA';(groups[k]=groups[k]||[]).push(s);});
  var order=R.lineages.slice(); if(groups['NA'])order.push('NA');
  var tot=R.samples.length||1;
  var comp='<div class="lincomp-bar">'+order.map(function(k){var n=(groups[k]||[]).length;if(!n)return '';var w=100*n/tot;
    return '<div class="lseg" style="width:'+w.toFixed(3)+'%;background:'+linColor(k=='NA'?null:k)+'" title="'+esc(k)+': '+n+' ('+w.toFixed(1)+'%)"></div>';}).join('')+'</div>'+
    '<div class="lincomp-lab">'+order.map(function(k){var n=(groups[k]||[]).length;if(!n)return '';
      return '<span class="lchip"><i style="background:'+linColor(k=='NA'?null:k)+'"></i>'+esc(k)+' <b>'+n+'</b></span>';}).join('')+'</div>';
  el('lincomp').innerHTML=comp;
  var head='<tr><th class="s" style="text-align:left">Lineage</th><th>n</th><th>%PASS</th><th>med depth</th><th>med breadth</th><th>med SNPs</th><th># MIXED</th></tr>';
  var body=order.map(function(k){var g=groups[k]||[];if(!g.length)return '';
    var np=g.filter(function(s){return s.v=='PASS';}).length,pct=g.length?100*np/g.length:0;
    var md=med(g.map(function(s){return s.m.mean_depth;}).filter(function(v){return v!=null;}));
    var mb=med(g.map(function(s){return s.m.breadth_pct;}).filter(function(v){return v!=null;}));
    var ms=med(g.map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}));
    var nmix=g.filter(function(s){return s.f.indexOf('MIXED')>=0;}).length;
    return '<tr data-lin="'+esc(k)+'"><td class="s" style="text-align:left"><span class="ldot" style="background:'+linColor(k=='NA'?null:k)+'"></span>'+esc(k)+'</td>'+
      '<td>'+g.length+'</td><td>'+pct.toFixed(0)+'</td>'+
      '<td>'+(md==null?'<span class="na">NA</span>':md.toFixed(1))+'</td><td>'+(mb==null?'<span class="na">NA</span>':mb.toFixed(1))+'</td>'+
      '<td>'+(ms==null?'<span class="na">NA</span>':Math.round(ms).toLocaleString('en-US'))+'</td>'+
      '<td'+(nmix>0?' style="color:var(--fail);font-weight:600"':'')+'>'+nmix+'</td></tr>';}).join('');
  el('linsumtable').innerHTML='<thead>'+head+'</thead><tbody>'+body+'</tbody>';
  Array.prototype.forEach.call(el('linsumtable').querySelectorAll('tbody tr'),function(tr){tr.onclick=function(){
    var k=tr.getAttribute('data-lin'); st.linFilter=(st.linFilter==(k=='NA'?null:k))?null:(k=='NA'?null:k);
    renderAll();el('gstats').scrollIntoView();};});
}

// per-sample genome sparkline for the detail modal: missing fraction (area) + SNP density (line) along the reference
function genomeSpark(s){
  var miss=s.miss, snp=s.trk&&s.trk.snp;
  if(!miss&&!snp)return '';
  var nb=(miss?miss.length:snp.length), gl=R.genome_len||nb;
  var W=512,H=46,pl=4,pr=4,pt=6,pb=12,iw=W-pl-pr,ih=H-pt-pb;
  function x(i){return pl+i/nb*iw;}
  var svg='<svg width="'+W+'" height="'+H+'" viewBox="0 0 '+W+' '+H+'" preserveAspectRatio="none" style="width:100%;height:'+H+'px">';
  svg+='<line x1="'+pl+'" y1="'+(pt+ih)+'" x2="'+(W-pr)+'" y2="'+(pt+ih)+'" stroke="'+TH.grid+'"/>';
  if(R.mask_bins)for(var mk=0;mk<nb&&mk<R.mask_bins.length;mk++){var mf=R.mask_bins[mk]||0;if(mf<0.5)continue;
    svg+='<rect x="'+x(mk).toFixed(1)+'" y="'+pt+'" width="'+(iw/nb+0.6).toFixed(1)+'" height="'+ih+'" fill="'+TH.faint+'" opacity="0.10"/>';}
  if(miss){var pts=[];for(var i=0;i<nb;i++){var mv=miss[i]==null?0:miss[i];pts.push(x(i).toFixed(1)+','+(pt+ih-mv*ih).toFixed(1));}
    svg+='<path d="M '+pl+' '+(pt+ih)+' L '+pts.join(' L ')+' L '+(W-pr)+' '+(pt+ih)+' Z" fill="rgba(214,64,58,.16)"/>'+
         '<polyline points="'+pts.join(' ')+'" fill="none" stroke="#d6403a" stroke-width="1"/>';}
  if(snp){var smax=Math.max.apply(null,snp.concat([1]));var sp=[];for(var j=0;j<nb;j++){sp.push(x(j).toFixed(1)+','+(pt+ih-(snp[j]||0)/smax*ih).toFixed(1));}
    svg+='<polyline points="'+sp.join(' ')+'" fill="none" stroke="#0e8ba8" stroke-width="1" opacity="0.75"/>';}
  [0,0.5,1].forEach(function(t){svg+='<text x="'+x(nb*t).toFixed(1)+'" y="'+(H-2)+'" font-size="7.5" fill="#94a3b8" text-anchor="'+(t==0?'start':t==1?'end':'middle')+'">'+fmtpos(Math.round(t*gl))+'</text>';});
  svg+='<text x="'+(W-pr)+'" y="'+(pt+6)+'" font-size="7.5" fill="#94a3b8" text-anchor="end">'+(snp?'missing (red) / SNP density (blue)':'missing along reference')+'</text>';
  return '<div class="dsub">Genome profile <span style="font-weight:400;color:#94a3b8">(is the loss one block or scattered?)</span></div><div style="margin-bottom:6px">'+svg+'</svg></div>';
}
// ---- aDNA damage authentication panel ----
function ctSpark(prof,lim){
  var W=132,H=34,pl=3,pr=3,pt=4,pb=8,iw=W-pl-pr,ih=H-pt-pb;
  if(!prof||!prof.length)return '<span class="nd" style="padding:0">no C&gt;T profile</span>';
  var n=Math.min(prof.length,25),xs=function(i){return pl+(n<=1?0:i/(n-1)*iw);};
  var ymax=Math.max(0.30,Math.max.apply(null,prof.slice(0,n).map(function(p){return p[1];})));
  var ys=function(v){return pt+ih-(v/ymax)*ih;};
  var pts=prof.slice(0,n).map(function(p,i){return xs(i).toFixed(1)+','+ys(p[1]).toFixed(1);}).join(' ');
  var tl=(lim!=null)?'<line x1="'+pl+'" y1="'+ys(lim).toFixed(1)+'" x2="'+(W-pr)+'" y2="'+ys(lim).toFixed(1)+'" stroke="#e0544f" stroke-width="1" stroke-dasharray="2 2"/>':'';
  return '<svg width="'+W+'" height="'+H+'" style="overflow:visible">'+
    '<line x1="'+pl+'" y1="'+(pt+ih)+'" x2="'+(W-pr)+'" y2="'+(pt+ih)+'" stroke="'+TH.grid+'"/>'+tl+
    '<polyline points="'+pts+'" fill="none" stroke="var(--accent)" stroke-width="1.5" stroke-linejoin="round"/>'+
    '<circle cx="'+xs(0).toFixed(1)+'" cy="'+ys(prof[0][1]).toFixed(1)+'" r="2.2" fill="var(--accent)"/>'+
    '<text x="'+pl+'" y="'+H+'" font-size="7.5" fill="#94a3b8">5&#39; pos</text>'+
    '<text x="'+(W-pr)+'" y="'+(pt-1)+'" font-size="7.5" fill="#94a3b8" text-anchor="end">'+(ymax*100).toFixed(0)+'%</text></svg>';
}
// ---- Functional gene burden (cohort HIGH+MODERATE snpEff burden per gene; optional --gene-burden aux TSV) ----
var GENEBURDEN_CAPTION='Genes ranked by cohort HIGH + MODERATE snpEff burden (summed across all samples); dominant effect = the most frequent effect class in that gene. HIGH + MODERATE counts are absolute, so a longer or more divergent gene ranks higher; this is a functional companion to the positional Variable genes panel below, not a selection test.';
function renderGeneBurden(){
  var host=el('gb_body'),tbl=el('gbtable'),note=el('gb_note'); if(!host||!tbl)return;
  var gb=R.gene_burden; if(!(gb&&gb.length)){host.style.display='none';return;} host.style.display='';
  var q=(st.gbq||'').toLowerCase(),exp=(el('geneburdenPanel')&&el('geneburdenPanel').classList.contains('expanded'));
  var flt=q?gb.filter(function(g){return (g.gene||'').toLowerCase().indexOf(q)>=0;}):gb;
  var rows=flt.slice(0,q?300:(exp?150:20));
  tbl.innerHTML='<thead><tr><th class="s" style="text-align:left">Gene</th><th>HIGH</th><th>MODERATE</th><th style="text-align:left">Dominant effect</th><th>Samples</th><th>Burden</th></tr></thead><tbody>'+
    rows.map(function(g){var mark=(g.start!=null&&g.end!=null);var burd=(g.total_impactful!=null?g.total_impactful:(g.high||0)+(g.moderate||0));
      return '<tr'+(mark?' data-start="'+g.start+'" data-end="'+g.end+'" data-name="'+esc(g.gene)+'" style="cursor:pointer"':'')+'>'+
        '<td class="s" style="text-align:left">'+esc(g.gene)+geneRvTag(g.gene)+'</td>'+
        '<td'+((g.high||0)>0?' style="color:var(--fail);font-weight:600"':'')+'>'+(g.high||0)+'</td>'+
        '<td>'+(g.moderate||0)+'</td>'+
        '<td style="text-align:left">'+esc(g.dominant_effect||'NA')+'</td>'+
        '<td>'+(g.n_samples==null?'-':g.n_samples)+'</td>'+
        '<td>'+Math.round(burd).toLocaleString('en-US')+'</td></tr>';}).join('')+'</tbody>';
  var gbcount=q?('Showing '+rows.length+' of '+flt.length+' matching genes. '):(flt.length>rows.length?('Top '+rows.length+' of '+flt.length+' genes (search or expand for more). '):'');
  note.textContent=gbcount+GENEBURDEN_CAPTION+(rows.some(function(g){return g.start!=null;})?'':' Click-to-mark is available when gene coordinates are joined from a GFF upstream.');
  if(!rows.length)tbl.innerHTML='<tbody><tr><td class="na" colspan="6" style="text-align:left;padding:8px">no gene matches "'+esc(st.gbq||'')+'".</td></tr></tbody>';
  Array.prototype.forEach.call(tbl.querySelectorAll('tbody tr[data-start]'),function(tr){tr.onclick=function(){
    var gl=R.genome_len||0,nb=R.nbins||200; if(!gl)return;
    var b0=Math.max(0,Math.min(nb-1,Math.floor(+tr.getAttribute('data-start')/gl*nb)));
    var b1=Math.max(0,Math.min(nb-1,Math.floor(+tr.getAttribute('data-end')/gl*nb)));
    st.gtrack='snp'; st.geneMark={b0:b0,b1:b1,name:tr.getAttribute('data-name')||''};
    renderGenome(); var gs=el('genome');if(gs)gs.scrollIntoView();};});
}
// ---- Selection: cohort per-gene dN/dS from eskaks (optional --pnps; population genetics, NOT per-sample QC) ----
var PNPS_CAPTION='Cohort selection analysis (pairwise dN/dS across samples per gene, eskaks Nei/Li). This is population genetics, NOT per-sample QC: a gene here is not a flag on any one sample. dN/dS > 1 is a signal of positive / diversifying selection at cohort scale but is noisy per gene, sensitive to alignment and to the codon model, and saturates at low divergence. dN, dS and the ratio can be undefined at low divergence (shown as NA when dS is near zero). This panel appears only when a cohort pN/pS table is provided upstream; it is not derivable from one annotated VCF.';
function renderPnps(){
  var host=el('pnps_body'),tbl=el('pnpstable'),note=el('pnps_note'); if(!host||!tbl)return;
  var pp=R.pnps; if(!(pp&&pp.length)){host.style.display='none';return;} host.style.display='';
  var q=(st.pnpsq||'').toLowerCase(),exp=(el('pnpsPanel')&&el('pnpsPanel').classList.contains('expanded'));
  var flt=q?pp.filter(function(g){return (g.gene||'').toLowerCase().indexOf(q)>=0;}):pp;
  var rows=flt.slice(0,q?300:(exp?150:40));
  tbl.innerHTML='<thead><tr><th class="s" style="text-align:left">Gene</th><th>dN</th><th>dS</th><th>dN/dS</th><th>pairs</th><th style="text-align:left">Effect</th></tr></thead><tbody>'+
    rows.map(function(g){var pos=(g.pnps!=null&&g.pnps>1);
      return '<tr><td class="s" style="text-align:left">'+esc(g.gene)+geneRvTag(g.gene)+'</td>'+
        '<td>'+(g.pn==null?'NA':g.pn.toFixed(3))+'</td>'+
        '<td>'+(g.ps==null?'NA':g.ps.toFixed(3))+'</td>'+
        '<td'+(pos?' style="color:var(--warn);font-weight:600"':'')+'>'+(g.pnps==null?'NA':g.pnps.toFixed(2))+'</td>'+
        '<td>'+(g.n==null?'-':Math.round(g.n))+'</td>'+
        '<td style="text-align:left">'+esc(g.effect||'NA')+'</td></tr>';}).join('')+'</tbody>';
  if(note)note.textContent=(q?('Showing '+rows.length+' of '+flt.length+' matching genes. '):'')+PNPS_CAPTION;
  if(!rows.length&&tbl)tbl.innerHTML='<tbody><tr><td class="na" colspan="6" style="text-align:left;padding:8px">no gene matches "'+esc(st.pnpsq||'')+'".</td></tr></tbody>';
}
function renderADNA(){
  var host=el('adna_body'); if(!host)return;
  var anc=R.samples.filter(function(s){return s.anc;});
  if(!R.n_ancient||!anc.length){host.innerHTML='<span class="nd" style="padding:0">no ancient (aDNA) samples in this run.</span>';return;}
  var lim=(athr&&athr.damage_min_ct!=null)?athr.damage_min_ct:null;
  var pctv=function(v){return v==null?'NA':(v*100).toFixed(1)+'%';};
  var cards=anc.map(function(s){var d=s.dmg,low=s.f.indexOf('DAMAGE_LOW')>=0;
    var ctcls=(d&&d.ct1!=null)?(low?'bad':'good'):'na';
    var head='<div class="acard-h"><span class="sname" data-s="'+esc(s.s)+'">'+esc(s.s)+'</span><span class="v '+s.v+'" style="font-size:10px">'+s.v+'</span></div>';
    if(!d)return '<div class="acard'+(low?' low':'')+'">'+head+'<div class="nd" style="padding:6px 0">damage data not found (run mapDamage2 / pass --mapdamage-dir)</div></div>';
    var body='<div class="acard-row"><div class="aspark">'+ctSpark(d.ct_profile,lim)+'</div><div class="astats">'+
      '<div class="astat"><span class="ak">5&#39; C&gt;T</span><span class="av '+ctcls+'">'+pctv(d.ct1)+'</span></div>'+
      '<div class="astat"><span class="ak">3&#39; G&gt;A</span><span class="av">'+pctv(d.ga1)+'</span></div>'+
      '<div class="astat"><span class="ak">frag len</span><span class="av">'+(d.fraglen==null?'NA':d.fraglen.toFixed(0)+' bp')+'</span></div></div></div>'+
      (low?'<div class="alow">'+icon('alert','sort')+'terminal C&gt;T below the '+(lim*100).toFixed(0)+'% authentication floor - possible modern contamination; verify before using as a calibration tip.</div>':'');
    return '<div class="acard'+(low?' low':'')+'">'+head+body+'</div>';}).join('');
  host.innerHTML='<div class="anote">Elevated terminal C&gt;T (5&#39;) / G&gt;A (3&#39;) deamination is <b>consistent with</b> post-mortem damage and is a necessary authentication signal. It does <b>not</b> by itself prove the DNA is ancient (deaminated modern contaminant DNA can mimic it) or endogenous, and its absence (below the floor) is a red flag for a modern sample mislabelled ancient. Treat this as a screen, not a proof.</div>'+
    '<div class="acards">'+cards+'</div>';
  Array.prototype.forEach.call(host.querySelectorAll('.sname'),function(sp){sp.onclick=function(){openDetail(sp.getAttribute('data-s'));};});
}

