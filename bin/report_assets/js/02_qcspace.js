var PCA_CAPTION='QC-metric space. Samples ordinated by their quality-metric profiles (each metric robustly z-scored across the samples in view, then projected onto the top two principal components of the metric covariance). This is a descriptive map of how QC profiles co-vary, not a map of biological relatedness and not a test. PC1 usually captures a coverage-and-completeness composite because depth, breadth and callability are strongly correlated; see the loadings for what drives each axis. A sample sitting apart can reflect genuine sequence divergence, a library or coverage artifact, or contamination, and the ordination alone cannot separate these, so cross-read the divergence-vs-completeness panel below. Axes are cohort-relative and recompute when you filter. The shaded ellipse marks the data spread (2 SD along each PC), not a confidence or significance region. Rings mark robust-Mahalanobis outliers.';
var REFBIAS_CAPTION='Reference divergence vs completeness. Each sample is plotted by consensus incompleteness (x) against a proxy for how far it sits from the reference (y = SNPs, per callable Mb when available to avoid coupling the two axes). This is a screen that separates two ways a sample can look reference-like. Bottom-left (few SNPs, little missing) is the reference-bias suspect zone: low apparent divergence that is not explained by missing data, consistent with reference-guided calling collapsing a divergent sample onto the reference, but a genuinely reference-close or smaller-genome sample looks identical here, so this is a list to check, not a verdict. Bottom-right (few SNPs, lots missing) is honest low coverage: the low count is just missing data, not bias. High divergence is not automatically real either, since contamination or a wrong reference also inflates SNP counts, so cross-read the mapping and error-rate metrics. Guide lines are the live missing gate (x) and a robust low-divergence cut (y = cohort median minus 2.5 robust SD); only samples below that cut are called out, so a normally-divergent sample is never flagged here. Confirm any suspected reference bias downstream; this panel only points.';
var TEMPORAL_CAPTION='Temporal sampling overview. Per group, the span and count of samples carrying a parseable collection year, alongside a rough sequence-information proxy (median SNPs vs reference). {DATED} of {TOT} samples carry a parseable year; overall span {SPAN}. This describes whether the cohort has the dated spread a downstream molecular-dating analysis would need, and it does not test for clock signal. Wide sampling in time is necessary but not sufficient: a cohort can span decades and still show no root-to-tip temporal signal, which only a proper date-randomization or regression test downstream can establish. The SNP median is an information proxy, not a count of clock-informative sites.';
// core quality/completeness metrics for the ordination + Mahalanobis. Divergence (snps/snp_density) is DELIBERATELY
// excluded here (a divergent lineage is not a QC outlier); it is the y-axis of the divergence-vs-completeness panel.
var PCA_CORE=['mean_depth','breadth_pct','callable_pct','missing_pct','iupac_pct','mapped_pct','duplication_pct','coverage_cv','mapq0_pct','error_rate_pct','q30_pct','het_frac','ti_tv','insert_size'];
function pcaColUsable(pool,k){var v=[],i;for(i=0;i<pool.length;i++){var x=pool[i].m[k];if(x!=null&&isFinite(x))v.push(x);}
  if(v.length<Math.max(4,Math.ceil(0.5*pool.length)))return null;   // require >50% present AND >=4 observations
  v.sort(function(a,b){return a-b;});var md=med2(v);
  var d=v.map(function(x){return Math.abs(x-md);}).sort(function(a,b){return a-b;});var sig=1.4826*med2(d);
  if(!(sig>0)){var mn=v.reduce(function(a,b){return a+b;},0)/v.length;sig=Math.sqrt(v.reduce(function(a,b){return a+(b-mn)*(b-mn);},0)/v.length);}
  if(!(sig>0))return null;return {med:md,sig:sig};}                 // truly constant -> drop the column
function buildZ(pool){var keys=[],meds=[],sigs=[],c;
  for(c=0;c<PCA_CORE.length;c++){var k=PCA_CORE[c],u=pcaColUsable(pool,k);if(u){keys.push(k);meds.push(u.med);sigs.push(u.sig);}}
  var K=keys.length,Z=[],samples=[],i,j;
  for(i=0;i<pool.length;i++){var s=pool[i],row=new Array(K);
    for(j=0;j<K;j++){var x=s.m[keys[j]];if(x==null||!isFinite(x))x=meds[j];var z=(x-meds[j])/sigs[j];if(z>8)z=8;else if(z<-8)z=-8;row[j]=z;}
    Z.push(row);samples.push(s);}
  return {Z:Z,keys:keys,samples:samples,k:K,n:Z.length};}
function covMatrix(Z,k){var n=Z.length,C=[],i,j,r;
  for(i=0;i<k;i++){C[i]=new Array(k);for(j=0;j<k;j++)C[i][j]=0;}
  for(r=0;r<n;r++){var row=Z[r];for(i=0;i<k;i++){var zi=row[i];for(j=i;j<k;j++)C[i][j]+=zi*row[j];}}
  var d=(n>1)?(n-1):1;for(i=0;i<k;i++)for(j=i;j<k;j++){C[i][j]/=d;C[j][i]=C[i][j];}return C;}
function topEigs(C,k,howMany){
  function matVec(M,x){var y=new Array(k),i,j;for(i=0;i<k;i++){var s=0;for(j=0;j<k;j++)s+=M[i][j]*x[j];y[i]=s;}return y;}
  function nrm(x){var s=0,i;for(i=0;i<k;i++)s+=x[i]*x[i];return Math.sqrt(s);}
  var A=[],a;for(a=0;a<k;a++)A[a]=C[a].slice();var vecs=[],vals=[],e,i,j;
  for(e=0;e<howMany;e++){var x=new Array(k);for(i=0;i<k;i++)x[i]=(((i*2654435761)>>>0)%1000)/1000-0.5+1e-6*(e+1);
    var nx=nrm(x)||1;for(i=0;i<k;i++)x[i]/=nx;var lam=0,it;
    for(it=0;it<500;it++){var y=matVec(A,x),ny=nrm(y);if(ny<1e-14)break;for(i=0;i<k;i++)y[i]/=ny;
      var Ay=matVec(A,y),nl=0;for(i=0;i<k;i++)nl+=y[i]*Ay[i];var dot=0;for(i=0;i<k;i++)dot+=y[i]*x[i];if(dot<0)for(i=0;i<k;i++)y[i]=-y[i];
      x=y;if(Math.abs(nl-lam)<1e-9*(Math.abs(nl)+1e-12)){lam=nl;break;}lam=nl;}
    if(lam<0)lam=0;vecs.push(x);vals.push(lam);
    for(i=0;i<k;i++)for(j=0;j<k;j++)A[i][j]-=lam*x[i]*x[j];}     // Hotelling deflation
  return {vecs:vecs,vals:vals};}
function pca2(pool){var B=buildZ(pool);
  if(B.n<4)return {ok:false,reason:'need at least 4 samples with a shared core metric in view',B:B};
  if(B.k<2)return {ok:false,reason:'need at least 2 core QC metrics in view',B:B};
  var C=covMatrix(B.Z,B.k),E=topEigs(C,B.k,2),tot=0,i;for(i=0;i<B.k;i++)tot+=C[i][i];
  if(!(tot>0))return {ok:false,reason:'no QC variance (samples identical on the core metrics)',B:B};
  var pev=[100*E.vals[0]/tot,100*E.vals[1]/tot];
  var scores=B.Z.map(function(row){var p1=0,p2=0,j;for(j=0;j<B.k;j++){p1+=row[j]*E.vecs[0][j];p2+=row[j]*E.vecs[1][j];}return [p1,p2];});
  function loadTop(vec){return B.keys.map(function(k,j){return {k:k,label:(MET[k]||{}).label||k,w:vec[j]};}).sort(function(a,b){return Math.abs(b.w)-Math.abs(a.w);}).slice(0,5);}
  return {ok:true,B:B,scores:scores,pev:pev,vecs:E.vecs,load:[loadTop(E.vecs[0]),loadTop(E.vecs[1])]};}
function invRidge(C,k,lambda){var A=[],I=[],i,j,r;
  for(i=0;i<k;i++){A[i]=C[i].slice();A[i][i]+=lambda;I[i]=new Array(k);for(j=0;j<k;j++)I[i][j]=(i==j)?1:0;}
  for(i=0;i<k;i++){var p=i,mx=Math.abs(A[i][i]);for(r=i+1;r<k;r++){var av=Math.abs(A[r][i]);if(av>mx){mx=av;p=r;}}
    if(mx<1e-12)return null;if(p!=i){var t=A[p];A[p]=A[i];A[i]=t;t=I[p];I[p]=I[i];I[i]=t;}
    var piv=A[i][i];for(j=0;j<k;j++){A[i][j]/=piv;I[i][j]/=piv;}
    for(r=0;r<k;r++){if(r==i)continue;var f=A[r][i];if(f===0)continue;for(j=0;j<k;j++){A[r][j]-=f*A[i][j];I[r][j]-=f*I[i][j];}}}
  return I;}
function robustCovInv(C,k){var tr=0,i;for(i=0;i<k;i++)tr+=C[i][i];var base=(tr/k)||1;var fac=[1e-3,1e-2,1e-1,1,10],t;
  for(t=0;t<fac.length;t++){var Iv=invRidge(C,k,fac[t]*base);if(Iv)return {inv:Iv,lambda:fac[t]*base};}return null;}
function chi2q975(k){var zp=1.959964,a=2/(9*k),tt=1-a+zp*Math.sqrt(a);return k*tt*tt*tt;}   // Wilson-Hilferty; nominal reference only
function mahalanobis(pool){var B=buildZ(pool);if(B.k<2||B.n<3)return {ok:false,B:B};
  var C=covMatrix(B.Z,B.k),R2=robustCovInv(C,B.k);if(!R2)return {ok:false,B:B};var Iv=R2.inv,rows=[],r;
  for(r=0;r<B.n;r++){var z=B.Z[r],acc=0,i,j;for(i=0;i<B.k;i++){var Iz=0;for(j=0;j<B.k;j++)Iz+=Iv[i][j]*z[j];acc+=z[i]*Iz;}if(acc<0)acc=0;
    var drv=B.keys.map(function(kk,ki){return {k:kk,label:(MET[kk]||{}).label||kk,z:z[ki]};}).sort(function(a,b){return Math.abs(b.z)-Math.abs(a.z);}).slice(0,3);
    rows.push({s:B.samples[r],d2:acc,drivers:drv});}
  return {ok:true,B:B,rows:rows};}
function spreadEllipse(pts){if(pts.length<4)return null;var n=pts.length,mx=0,my=0,i;   // 2-SD spread ellipse in pixel space
  for(i=0;i<n;i++){mx+=pts[i][0];my+=pts[i][1];}mx/=n;my/=n;
  var sxx=0,syy=0,sxy=0;for(i=0;i<n;i++){var dx=pts[i][0]-mx,dy=pts[i][1]-my;sxx+=dx*dx;syy+=dy*dy;sxy+=dx*dy;}
  var d=(n-1)||1;sxx/=d;syy/=d;sxy/=d;var tr=sxx+syy,det=sxx*syy-sxy*sxy,disc=Math.sqrt(Math.max(0,tr*tr/4-det));
  var l1=tr/2+disc,l2=tr/2-disc;if(l1<0)l1=0;if(l2<0)l2=0;
  var ang=(Math.abs(sxy)<1e-12&&sxx>=syy)?0:Math.atan2(l1-sxx,sxy);
  return {cx:mx,cy:my,rx:2*Math.sqrt(l1),ry:2*Math.sqrt(l2),angle:ang};}
function yearOf(dstr){if(dstr==null)return null;var m=String(dstr).match(/(1[5-9]\d\d|20\d\d|2100)/);if(!m)return null;var y=+m[1];return (y>=1500&&y<=2100)?y:null;}
// register d^2 as a first-class metric + a WARN-only flag (never in FAILF -> never a FAIL, never auto-excluded)
R.metrics.push({key:'qc_mahal',label:'QC outlier d²',kind:'float',dir:'hi_bad'});
if(DIST.indexOf('qc_mahal')<0)DIST.push('qc_mahal');
R.defs.qc_mahal=['Robust Mahalanobis distance squared in standardized QC-metric space (ridge-regularized robust covariance, median center). A multivariate rarity rank: high = an unusual COMBINATION of QC metrics even when each metric is individually in range. A heuristic composite, not a hypothesis test; missing metrics are imputed to the cohort median so a sparse sample reads more central than it may be.','robust flag QC_OUTLIER when d2 > median + 3*MAD of the cohort d2 distribution'];
R.mahal_cut=null;R.mahal_k=null;R.mahal_ref=null;R.mahal_rows=null;
// compute d^2 over the WHOLE cohort (stable flag under filtering); called inside recompute() so it stays live.
function applyMahal(){
  var pool=R.samples.filter(function(s){return !s.anc&&(!s.f||s.f.indexOf('NO_DATA')<0);});   // modern, produced-something (aDNA is legitimately extreme; scope like the SNP-z)
  R.samples.forEach(function(s){s.m.qc_mahal=null;s._mdrv=null;});                             // default: not scored (aDNA / NO_DATA -> NA, no ring, no flag)
  var M=mahalanobis(pool);
  if(!M.ok){R.mahal_cut=null;R.mahal_k=null;R.mahal_ref=null;R.mahal_rows=null;return;}
  var d2s=M.rows.map(function(o){return o.d2;}).sort(function(a,b){return a-b;});
  var md=med2(d2s),mad=med2(d2s.map(function(v){return Math.abs(v-md);}).sort(function(a,b){return a-b;}));
  var scale=1.4826*mad;   // MAD=0 (ties from median-imputation) -> fall back to SD; genuinely zero spread -> Infinity (flag nothing), like the SNP-z guard
  if(!(scale>0)){var mn=d2s.reduce(function(a,b){return a+b;},0)/d2s.length;scale=Math.sqrt(d2s.reduce(function(a,b){return a+(b-mn)*(b-mn);},0)/d2s.length);}
  R.mahal_cut=(scale>0)?md+3*scale:Infinity;R.mahal_k=M.B.k;R.mahal_ref=chi2q975(M.B.k);R.mahal_rows=M.rows;
  var by={};M.rows.forEach(function(o){by[o.s.s]=o;});
  R.samples.forEach(function(s){var o=by[s.s];if(o){s.m.qc_mahal=o.d2;s._mdrv=o.drivers;}});
  // refresh this live-derived metric's display scaling (RANGES/MED/IQR are built once at load, before d2 exists)
  var vs=R.samples.map(function(s){return s.m.qc_mahal;}).filter(function(v){return v!=null;});
  if(vs.length){RANGES.qc_mahal=[Math.min.apply(null,vs),Math.max.apply(null,vs)];var so=vs.slice().sort(function(a,b){return a-b;});MED.qc_mahal=so[Math.floor(so.length/2)];IQR.qc_mahal=[quant(so,0.25),quant(so,0.5),quant(so,0.75)];}}
// ---- Panel 1: QC-metric-space PCA biplot (over the live cohort) ----
function renderQCspace(){var host=el('qcpca_body'),cap=el('qcpca_caption'),ot=el('pca_outtable');if(!host)return;
  var P=pca2(visible());
  if(!P.ok){host.innerHTML='<div class="nd">'+esc(P.reason)+'</div>';if(cap)cap.innerHTML='';if(ot)ot.innerHTML='';return;}
  var pad=44,W=Math.min(host.clientWidth||620,760),H=Math.min(W,420);
  var xs=P.scores.map(function(p){return p[0];}),ys=P.scores.map(function(p){return p[1];});
  var xr=[Math.min.apply(null,xs),Math.max.apply(null,xs)],yr=[Math.min.apply(null,ys),Math.max.apply(null,ys)];
  var mgx=(xr[1]-xr[0])*0.08||1,mgy=(yr[1]-yr[0])*0.08||1;xr=[xr[0]-mgx,xr[1]+mgx];yr=[yr[0]-mgy,yr[1]+mgy];
  function sx(v){return pad+(v-xr[0])/(xr[1]-xr[0]||1)*(W-pad-14);}
  function sy(v){return H-pad-(v-yr[0])/(yr[1]-yr[0]||1)*(H-pad-16);}
  var groups={};P.B.samples.forEach(function(s,i){var g=(st.colorBy=='lineage')?(s.lineage||'NA'):s.v;(groups[g]=groups[g]||[]).push([sx(P.scores[i][0]),sy(P.scores[i][1])]);});
  var ellSVG='',gn;for(gn in groups){var E=spreadEllipse(groups[gn]);if(!E)continue;var col=(st.colorBy=='lineage')?linColor(gn):(VCOL[gn]||'#94a3b8');
    ellSVG+='<ellipse cx="'+E.cx.toFixed(1)+'" cy="'+E.cy.toFixed(1)+'" rx="'+E.rx.toFixed(1)+'" ry="'+E.ry.toFixed(1)+'" transform="rotate('+(E.angle*180/Math.PI).toFixed(1)+' '+E.cx.toFixed(1)+' '+E.cy.toFixed(1)+')" fill="'+col+'" fill-opacity="0.06" stroke="'+col+'" stroke-opacity="0.4" stroke-width="1"/>';}
  var dots=P.B.samples.map(function(s,i){var big=(st.hi==s.s),cx=sx(P.scores[i][0]).toFixed(1),cy=sy(P.scores[i][1]).toFixed(1);
    var ring=(R.mahal_cut!=null&&s.m.qc_mahal!=null&&s.m.qc_mahal>R.mahal_cut)?'<circle cx="'+cx+'" cy="'+cy+'" r="7" fill="none" stroke="'+VCOL.WARN+'" stroke-width="1.3"/>':'';
    return ring+'<circle cx="'+cx+'" cy="'+cy+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="0.85"'+((s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+P.scores[i][0]+'" data-y="'+P.scores[i][1]+'" data-xl="PC1" data-yl="PC2" data-xk="float" data-yk="float"/>';}).join('');
  var frame='<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(W-14)+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/><line x1="'+pad+'" y1="14" x2="'+pad+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/>';
  var zero='';if(xr[0]<0&&xr[1]>0)zero+='<line x1="'+sx(0).toFixed(1)+'" y1="14" x2="'+sx(0).toFixed(1)+'" y2="'+(H-pad)+'" stroke="#eef2f6"/>';if(yr[0]<0&&yr[1]>0)zero+='<line x1="'+pad+'" y1="'+sy(0).toFixed(1)+'" x2="'+(W-14)+'" y2="'+sy(0).toFixed(1)+'" stroke="#eef2f6"/>';
  var xt='<text x="'+((pad+W-14)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">PC1 ('+P.pev[0].toFixed(1)+'%)</text>';
  var yt='<text transform="rotate(-90 13 '+((14+H-pad)/2)+')" x="13" y="'+((14+H-pad)/2)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">PC2 ('+(P.pev[1]<1e-3?'~0%, rank-deficient':P.pev[1].toFixed(1)+'%')+')</text>';
  var load='<div style="display:flex;gap:24px;flex-wrap:wrap;margin-top:8px">'+[0,1].map(function(pc){return '<div style="flex:1;min-width:170px"><div class="dsub" style="margin:2px 0 6px">PC'+(pc+1)+' loadings</div>'+P.load[pc].map(function(l){var w=Math.abs(l.w),col=l.w>=0?'#22a06b':'#e0544f';return '<div class="drow"><span class="dk">'+esc(l.label)+'</span><div class="dbarwrap"><div class="dbar" style="width:'+Math.round(w*100)+'%;background:'+col+'"></div></div><span class="dv">'+(l.w>=0?'+':'')+l.w.toFixed(2)+'</span></div>';}).join('')+'</div>';}).join('')+'</div>';
  host.innerHTML='<svg viewBox="0 0 '+W+' '+H+'" style="width:100%;height:auto;display:block;cursor:crosshair">'+ellSVG+frame+zero+xt+yt+dots+'</svg>'+colorLegend()+load;
  if(cap)cap.innerHTML=PCA_CAPTION;
  if(ot){var uv=visible().filter(function(s){return s.m.qc_mahal!=null;}).sort(function(a,b){return b.m.qc_mahal-a.m.qc_mahal;}).slice(0,8);
    ot.innerHTML='<thead><tr><th class="s" style="text-align:left">Sample</th><th>d²</th><th>QC</th><th style="text-align:left">Top deviating metrics</th></tr></thead><tbody>'+
     (uv.length?uv.map(function(s){var drv=(s._mdrv||[]).map(function(dd){var col=dd.z>=0?'#a01f2d':'#0b7350';return '<span style="color:'+col+'">'+esc(dd.label)+' '+(dd.z>=0?'+':'')+dd.z.toFixed(1)+'σ</span>';}).join(', ');
       return '<tr data-s="'+esc(s.s)+'"'+(st.hi==s.s?' class="hl"':'')+'><td class="s"><span class="sname" data-s="'+esc(s.s)+'">'+esc(s.s)+'</span></td><td>'+shortv(s.m.qc_mahal,'float')+'</td><td><span class="v '+s.v+'">'+s.v+'</span></td><td style="text-align:left;white-space:normal">'+drv+'</td></tr>';}).join('')
      :'<tr><td colspan="4" style="text-align:left;color:#94a3b8;padding:8px">no samples in view</td></tr>')+'</tbody>';
    Array.prototype.forEach.call(ot.querySelectorAll('tbody tr[data-s]'),function(tr){tr.onclick=function(){setHi(tr.getAttribute('data-s'));};});}}
// ---- Panel 2: divergence vs completeness (reference-bias screen) ----
function renderRefBias(){var host=el('divcomp_body'),cap=el('divcomp_caption'),qn=el('divcomp_quadn');if(!host)return;
  st.divx=st.divx||'missing_pct';
  var yk=(R.snp_density_ok?'snp_density':'snps'),ylabel=(R.snp_density_ok?'SNPs per callable Mb':'SNPs vs reference');
  function xof(s){return st.divx=='callable_inv'?(s.m.callable_pct!=null?100-s.m.callable_pct:null):(s.m.missing_pct!=null?s.m.missing_pct:(s.m.callable_pct!=null?100-s.m.callable_pct:null));}
  var pool=visible().filter(function(s){return xof(s)!=null&&s.m[yk]!=null;});
  if(pool.length<3){host.innerHTML='<div class="nd">needs missing/callable and SNP data for at least 3 samples in view</div>';if(cap)cap.innerHTML='';if(qn)qn.innerHTML='';return;}
  var pad=44,W=Math.min(host.clientWidth||620,760),H=Math.min(W*0.72,360);
  var xs=pool.map(xof),ys=pool.map(function(s){return s.m[yk];});
  var xhi=Math.max.apply(null,xs)*1.06||1,yhi=Math.max.apply(null,ys)*1.06||1;
  function sx(v){return pad+v/(xhi||1)*(W-pad-14);}function sy(v){return H-pad-v/(yhi||1)*(H-pad-16);}
  var gx=thr.missing_max,divs=ys.slice().sort(function(a,b){return a-b;});
  var dmed=med2(divs),dmad=med2(divs.map(function(v){return Math.abs(v-dmed);}).sort(function(a,b){return a-b;}));
  var ylo=dmed-2.5*1.4826*dmad;if(!(dmad>0)||!(ylo>0))ylo=0.5*dmed;   // robust low cut (guard MAD=0 ties + wide spread): only ANOMALOUSLY low divergence is called out (not just below median)
  var GX=sx(Math.min(gx,xhi)),GY=sy(Math.max(0,ylo));
  var rects='<rect x="'+pad+'" y="14" width="'+Math.max(0,W-14-pad).toFixed(1)+'" height="'+Math.max(0,GY-14).toFixed(1)+'" fill="var(--pass)" opacity="0.04"/>'
    +'<rect x="'+pad+'" y="'+GY.toFixed(1)+'" width="'+Math.max(0,GX-pad).toFixed(1)+'" height="'+Math.max(0,H-pad-GY).toFixed(1)+'" fill="var(--fail)" opacity="0.07"/>'
    +'<rect x="'+GX.toFixed(1)+'" y="'+GY.toFixed(1)+'" width="'+Math.max(0,W-14-GX).toFixed(1)+'" height="'+Math.max(0,H-pad-GY).toFixed(1)+'" fill="var(--warn)" opacity="0.06"/>';
  var guides='<line x1="'+GX.toFixed(1)+'" y1="14" x2="'+GX.toFixed(1)+'" y2="'+(H-pad)+'" stroke="#94a3b8" stroke-dasharray="4 3"/><line x1="'+pad+'" y1="'+GY.toFixed(1)+'" x2="'+(W-14)+'" y2="'+GY.toFixed(1)+'" stroke="#94a3b8" stroke-dasharray="4 3"/>';
  var quad={refbias:0,lowcov:0,typical:0};
  var dots=pool.map(function(s){var xv=xof(s),yv=s.m[yk],lowMiss=xv<=gx,lowDiv=yv<ylo;
    var cls=lowDiv?(lowMiss?'refbias':'lowcov'):'typical';quad[cls]++;
    var big=(st.hi==s.s),cx=sx(xv).toFixed(1),cy=sy(yv).toFixed(1);
    var ring=(cls=='refbias')?'<circle cx="'+cx+'" cy="'+cy+'" r="6.5" fill="none" stroke="#a01f2d" stroke-width="1.2"/>':'';
    return ring+'<circle cx="'+cx+'" cy="'+cy+'" r="'+(big?5.4:3.4)+'" fill="'+dotColor(s)+'" opacity="0.82"'+((s.v=='FAIL'||big)?' stroke="'+TH.ink+'" stroke-width="'+(big?1.4:0.6)+'"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+xv+'" data-y="'+yv+'" data-xl="'+(st.divx=='callable_inv'?'100 - callable %':'Missing %')+'" data-yl="'+esc(ylabel)+'" data-xk="pct" data-yk="'+(R.snp_density_ok?'float':'int')+'"/>';}).join('');
  var frame='<line x1="'+pad+'" y1="'+(H-pad)+'" x2="'+(W-14)+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/><line x1="'+pad+'" y1="14" x2="'+pad+'" y2="'+(H-pad)+'" stroke="'+TH.axis+'"/>';
  var titles='<text x="'+((pad+W-14)/2)+'" y="'+(H-6)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">'+(st.divx=='callable_inv'?'100 - callable % (incompleteness)':'Missing % (incompleteness)')+'</text><text transform="rotate(-90 13 '+((14+H-pad)/2)+')" x="13" y="'+((14+H-pad)/2)+'" text-anchor="middle" font-size="11" fill="'+TH.mut+'">'+esc(ylabel)+'</text>';
  var labs='<text x="'+(pad+6)+'" y="'+(H-pad-6)+'" font-size="8.5" font-weight="600" fill="var(--fail)">reference-bias suspect</text>'
    +'<text x="'+(pad+6)+'" y="24" font-size="8.5" font-weight="600" fill="#3f7d55">typical divergence</text>'
    +'<text x="'+(W-16)+'" y="'+(H-pad-6)+'" text-anchor="end" font-size="8.5" font-weight="600" fill="var(--warn)">low coverage</text>';
  host.innerHTML='<svg viewBox="0 0 '+W+' '+H+'" style="width:100%;height:auto;display:block;cursor:crosshair">'+rects+guides+frame+titles+labs+dots+'</svg>'+colorLegend();
  if(qn)qn.innerHTML=[['refbias','reference-bias suspect','var(--fail)'],['lowcov','honest low-coverage','var(--warn)'],['typical','typical divergence','var(--pass)']].map(function(q){return '<span><i style="background:'+q[2]+'"></i>'+q[1]+' <b>'+quad[q[0]]+'</b></span>';}).join('');
  if(cap)cap.innerHTML=REFBIAS_CAPTION;
  var dxb=el('divx');if(dxb)Array.prototype.forEach.call(dxb.querySelectorAll('button'),function(b){b.onclick=function(){st.divx=b.getAttribute('data-x');Array.prototype.forEach.call(dxb.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});renderRefBias();};});}
// ---- Panel 3: temporal sampling overview (per lineage) ----
function renderTemporal(){var host=el('temporal_body'),cap=el('temporal_caption');if(!host)return;
  var pool=visible();
  var dated=pool.map(function(s){return {s:s,y:yearOf(s.date)};}).filter(function(o){return o.y!=null;});
  if(!dated.length){host.innerHTML='<div class="anote">Temporal overview needs sample collection years (a <b>date</b> column); none parseable in this run. A time-resolved / clock analysis is a downstream step.</div>';if(cap)cap.innerHTML='';return;}
  var ymin=Math.min.apply(null,dated.map(function(o){return o.y;})),ymax=Math.max.apply(null,dated.map(function(o){return o.y;})),span=Math.max(1,ymax-ymin);
  var g={};pool.forEach(function(s){var L=R.lin_present?(s.lineage||'NA'):'all samples';(g[L]=g[L]||[]).push(s);});
  var order=R.lin_present?(R.lineages||[]).concat(['NA']):['all samples'],seen={},rows=[];
  order.forEach(function(L){if(g[L]&&!seen[L]){seen[L]=1;rows.push(L);}});Object.keys(g).forEach(function(L){if(!seen[L]){seen[L]=1;rows.push(L);}});
  var gut=140,W=Math.min(host.clientWidth||620,900),plotW=Math.max(80,W-gut-70);
  function ax(y){return gut+(y-ymin)/span*plotW;}
  var ticks='';for(var t=0;t<=4;t++){var yr=Math.round(ymin+span*t/4),X=ax(yr);ticks+='<line x1="'+X.toFixed(1)+'" y1="14" x2="'+X.toFixed(1)+'" y2="18" stroke="'+TH.axis+'"/><text x="'+X.toFixed(1)+'" y="11" text-anchor="middle" font-size="8.5" fill="#94a3b8">'+yr+'</text>';}
  var axisSVG='<svg viewBox="0 0 '+W+' 22" style="width:100%;height:auto;display:block"><line x1="'+gut+'" y1="18" x2="'+(gut+plotW)+'" y2="18" stroke="'+TH.grid+'"/>'+ticks+'</svg>';
  var body=rows.map(function(L){var grp=g[L],ys=grp.map(function(s){return yearOf(s.date);}).filter(function(y){return y!=null;});
    var snps=grp.map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
    var titv=grp.map(function(s){return s.m.ti_tv;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
    var mn=ys.length?Math.min.apply(null,ys):null,mxx=ys.length?Math.max.apply(null,ys):null;
    var col=R.lin_present?linColor(L):'var(--accent)';
    var rug=grp.map(function(s,i){var y=yearOf(s.date);if(y==null)return '';var jx=(((i*2654435761)>>>0)%997)/997-0.5,big=(st.hi==s.s);
      return '<circle cx="'+ax(y).toFixed(1)+'" cy="'+(13+jx*8).toFixed(1)+'" r="'+(big?4.5:2.4)+'" fill="'+col+'" opacity="0.72"'+(big?' stroke="'+TH.ink+'" stroke-width="1"':'')+' data-s="'+esc(s.s)+'" data-lin="'+esc(s.lineage||'')+'" data-x="'+y+'" data-y="'+(s.m.snps!=null?s.m.snps:0)+'" data-xl="year" data-yl="SNPs" data-xk="int" data-yk="int"/>';}).join('');
    var bar=(mn!=null)?'<rect x="'+ax(mn).toFixed(1)+'" y="10" width="'+Math.max(2,ax(mxx)-ax(mn)).toFixed(1)+'" height="6" rx="3" fill="'+col+'" opacity="0.18"/>':'';
    var lineSVG='<svg viewBox="0 0 '+W+' 26" style="width:100%;height:auto;display:block">'+bar+rug+'</svg>';
    var proxy=snps.length?med2(snps):null,mtitv=titv.length?med2(titv):null;
    return '<div class="bee" style="height:auto;min-height:34px"><div class="bl"><span class="ldot" style="background:'+col+'"></span>'+esc(L)+'<span class="stat">n='+grp.length+(mn!=null?' · '+mn+'-'+mxx:' · no dates')+(mtitv!=null?' · Ti/Tv '+mtitv.toFixed(2):'')+'</span></div><div class="plotwrap">'+lineSVG+'</div><div class="sv" style="flex:0 0 62px" title="informative-site proxy (median SNPs)">'+(proxy!=null?Math.round(proxy).toLocaleString('en-US'):'-')+'</div></div>';}).join('');
  host.innerHTML=axisSVG+body;
  var cohTitv=(function(){var tv=pool.map(function(s){return s.m.ti_tv;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});return tv.length?med2(tv):null;})();
  if(cap)cap.innerHTML=TEMPORAL_CAPTION.replace('{DATED}',dated.length).replace('{TOT}',pool.length).replace('{SPAN}',ymin+' - '+ymax+' ('+span+' yr)')+(cohTitv!=null?' Cohort median Ti/Tv = '+cohTitv.toFixed(2)+' (a ratio); a ratio drifting toward parity indicates mutational saturation where the SNP proxy over-counts, and a high ratio is spectrum or damage, not saturation.':'');}
// ===============================================================================
// ---- auto-discovered extra metrics: merge into the registry, hidden by default ----
