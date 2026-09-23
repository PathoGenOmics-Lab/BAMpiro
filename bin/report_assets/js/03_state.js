var extraSet={}; R.extra.forEach(function(e){extraSet[e.key]=1;R.metrics.push(e);R.defs[e.key]=R.defs[e.key]||['Auto-detected metric from the summary TSV (not a named QC metric).',''];});
var st={sortKey:'v',asc:false,q:'',onlyFlagged:false,hidden:{},hi:null,flagFilter:null,ptype:'beeswarm',   // default view: worst QC first (a returning user's saved sort overrides this)
        sx:'mean_depth',sy:'breadth_pct',excl:{},detail:null,colorBy:'qc',groupLin:false,ancOnly:null,linFilter:null,metaFilter:{},gtrack:'missing',maskOn:false,gsel:null,gzoom:null,gbq:'',hotq:'',pnpsq:'',vardoseq:'',colf:{},showColF:false};
var SGEO=null, GGEO=null;   // scatter + genome brush geometry caches (for inverse-mapping the rubber-band)
var _thdb,_qdb,_mxdb,_cfdb;  // debounce timers: keep live inputs snappy at cohort scale (defer heavy re-renders)
// The table opens on the metrics a verdict is made of; the other ~35 are one click away under "columns".
// Showing all of them at once made a 48-column table in which nothing could be read.
// They are listed in reading order: how much sequence, how much of the genome it covers, how usable the
// consensus is, then what was found in it; any other column the reader turns on follows them.
var TABLE_ORDER=['raw_reads','mean_depth','breadth_pct','breadth10x_pct','mapped_pct','duplication_pct','missing_pct',
  'callable_pct','iupac_pct','snps','het_frac','ti_tv','ann_high','n_lineages'];
var TABLE_KEY={}; TABLE_ORDER.forEach(function(k,i){TABLE_KEY[k]=i+1;});
R.metrics.forEach(function(m){if(!TABLE_KEY[m.key])st.hidden[m.key]=1;});
function tableMetrics(){return R.metrics.filter(function(m){return !st.hidden[m.key];}).map(function(m,i){return [m,i];})
  .sort(function(a,b){return ((TABLE_KEY[a[0].key]||99)-(TABLE_KEY[b[0].key]||99))||(a[1]-b[1]);}).map(function(x){return x[0];});}
// narrow (embedded panel / mobile): fewer still; the rest stays under "columns"
if(window.innerWidth<860){['raw_reads','breadth10x_pct','iupac_pct','het_frac','n_lineages','ann_high','callable_pct'].forEach(function(k){st.hidden[k]=1;});}
var MET={}; R.metrics.forEach(function(m){MET[m.key]=m;});
function quant(so,p){if(!so.length)return null;var i=(so.length-1)*p,lo=Math.floor(i),hi=Math.ceil(i);return lo==hi?so[lo]:so[lo]+(so[hi]-so[lo])*(i-lo);}
var RANGES={},MED={},IQR={}; R.metrics.forEach(function(m){
  var vs=R.samples.map(function(s){return s.m[m.key];}).filter(function(v){return v!=null;});
  RANGES[m.key]=vs.length?[Math.min.apply(null,vs),Math.max.apply(null,vs)]:[0,1];
  var so=vs.slice().sort(function(a,b){return a-b;}); MED[m.key]=so.length?so[Math.floor(so.length/2)]:null;
  IQR[m.key]=so.length?[quant(so,0.25),quant(so,0.5),quant(so,0.75)]:null;
});
// robust snp median/mad over MODERN samples (matches the Python cohort scoping)
// even-n midpoint average to match the Python robust() convention (so qc_flags.tsv and the HTML agree)
function med2(a){var n=a.length;return n%2?a[(n-1)/2]:(a[n/2-1]+a[n/2])/2;}
var SNPMED=null,SNPSIG=null;(function(){var vs=R.samples.filter(function(s){return !s.anc;}).map(function(s){return s.m.snps;}).filter(function(v){return v!=null;}).sort(function(a,b){return a-b;});
  if(!vs.length)return; var n=vs.length,med=med2(vs);var d=vs.map(function(v){return Math.abs(v-med);}).sort(function(a,b){return a-b;});
  var mad=med2(d),sig=1.4826*mad; if(sig<=0){var mean=vs.reduce(function(a,b){return a+b;},0)/n; sig=Math.sqrt(vs.reduce(function(a,b){return a+(b-mean)*(b-mean);},0)/n)||1;} SNPMED=med;SNPSIG=sig;})();
// ---- snpEff functional-annotation cohort fits (over NON-ANCIENT samples; robust median + MAD like the SNP-z) ----
function robustMedSig(vals){var v=vals.filter(function(x){return x!=null;}).sort(function(a,b){return a-b;});
  if(v.length<4)return null;var md=med2(v);
  var d=v.map(function(x){return Math.abs(x-md);}).sort(function(a,b){return a-b;});var sig=1.4826*med2(d);
  if(sig<=0){var mn=v.reduce(function(a,b){return a+b;},0)/v.length;sig=Math.sqrt(v.reduce(function(a,b){return a+(b-mn)*(b-mn);},0)/v.length)||0;}
  return sig>0?{med:md,sig:sig}:null;}
var ANNHI=null,LOFHI=null,ANNPCT=null,MSIL=null;
(function(){var mod=R.samples.filter(function(s){return !s.anc;});
  ANNHI =robustMedSig(mod.map(function(s){return s.m.ann_high;}));
  LOFHI =robustMedSig(mod.map(function(s){return s.m.lof_pct;}));
  ANNPCT=robustMedSig(mod.map(function(s){return s.m.annotated_pct;}));
  MSIL  =robustMedSig(mod.map(function(s){return s.m.missense_silent;}));})();
var ANNOT_FLOOR=(typeof thr.annot_floor=='number')?thr.annot_floor:50;   // absolute annotated-% sanity floor
var FAILF={LOW_DEPTH:1,LOW_BREADTH:1,HIGH_MISSING:1,NO_DATA:1};
// ---- expected-range bands: acceptable [lo,hi] per metric, derived LIVE from the ACTIVE threshold set ----
var TITV_WINDOW=[1.5,2.1];  // Ti/Tv sanity window (a ratio); A4/M. bovis legitimately runs higher (~2.6)
function bandFor(pk,T){switch(pk){
  case 'mean_depth':return [T.depth_min,Infinity];
  case 'breadth_pct':return [T.breadth_min,Infinity];
  case 'mapped_pct':return [T.mapping_min,Infinity];
  case 'missing_pct':return [-Infinity,T.missing_max];
  case 'duplication_pct':return [-Infinity,T.dup_max];
  case 'iupac_pct':return [-Infinity,T.iupac_max];
  case 'het_frac':return [-Infinity,thr.het_max_frac];
  case 'ti_tv':return [(thr.titv_min>0?Math.max(thr.titv_min,TITV_WINDOW[0]):TITV_WINDOW[0]),TITV_WINDOW[1]];
  case 'snps':return (SNPMED==null||!SNPSIG)?null:[SNPMED-thr.snp_z*SNPSIG,SNPMED+thr.snp_z*SNPSIG];
  case 'ann_high':return (ANNHI==null)?null:[-Infinity,ANNHI.med+thr.snp_z*ANNHI.sig];
  case 'annotated_pct':return (ANNPCT==null)?null:[ANNPCT.med-thr.snp_z*ANNPCT.sig,Infinity];
  case 'missense_silent':return null;   // deliberately no band: it is a proxy, high-tail only
  default:return null;}}
var RULES=[['LOW_DEPTH','mean_depth','<','depth_min','FAIL'],['LOW_BREADTH','breadth_pct','<','breadth_min','FAIL'],
 ['HIGH_MISSING','missing_pct','>','missing_max','FAIL'],['MAPPING_LOW','mapped_pct','<','mapping_min','WARN'],
 ['HIGH_DUP','duplication_pct','>','dup_max','WARN'],['HIGH_IUPAC','iupac_pct','>','iupac_max','WARN'],
 ['TITV_LOW','ti_tv','<','titv_min','WARN']];
// The flags recompute() derives itself from the live thresholds. Any other flag on a sample was
// decided by the pipeline from evidence the page does not hold (LINEAGE_MISMATCH compares the
// sample with the rest of its reference), so re-flagging keeps it instead of silently dropping it:
// dropping it hid every wrong-reference sample and let the page disagree with qc_flags.tsv.
var JS_FLAGS={LOW_DEPTH:1,LOW_BREADTH:1,HIGH_MISSING:1,NO_DATA:1,MAPPING_LOW:1,HIGH_DUP:1,HIGH_IUPAC:1,TITV_LOW:1,
  SNP_LOW:1,SNP_HIGH:1,HET_HIGH:1,HIGH_IMPACT_EXCESS:1,LOF_EXCESS:1,PNPS_PROXY_HIGH:1,ANNOTATION_POOR:1,MIXED:1,
  DAMAGE_LOW:1,QC_OUTLIER:1};
function pipelineFlags(s){return (s.f0||[]).filter(function(x){return !JS_FLAGS[x];});}
// What each flag means in words, for the chips and the flagged list (the code stays in the exports).
var FLAGLAB={LOW_DEPTH:'Low depth',LOW_BREADTH:'Low breadth',HIGH_MISSING:'Incomplete consensus',NO_DATA:'No output',
  MAPPING_LOW:'Few reads map',HIGH_DUP:'High duplication',HIGH_IUPAC:'Many ambiguous bases',TITV_LOW:'Low Ti/Tv',
  SNP_LOW:'Unusually few SNPs',SNP_HIGH:'Unusually many SNPs',HET_HIGH:'Mixed alleles',MIXED:'Mixed lineages',
  LINEAGE_MISMATCH:'Wrong reference?',DAMAGE_LOW:'No aDNA damage',QC_OUTLIER:'Unusual QC profile',
  HIGH_IMPACT_EXCESS:'Excess high-impact variants',LOF_EXCESS:'Excess loss-of-function',
  PNPS_PROXY_HIGH:'High missense/silent',ANNOTATION_POOR:'Poor annotation'};
function flagLab(f){return FLAGLAB[f]||String(f).toLowerCase().replace(/_/g,' ');}
function recompute(){
  R.samples.forEach(function(s){
    if(!s.f0)s.f0=(s.f||[]).slice();   // the pipeline's own verdict, kept for pipelineFlags()
    if(['mean_depth','breadth_pct','missing_pct','mapped_pct'].every(function(k){return s.m[k]==null;})){s.f=['NO_DATA'];s.v='FAIL';return;}
    var f=[],T=actv(s);
    RULES.forEach(function(r){var v=s.m[r[1]];if(v==null)return;
      var t=(r[0]=='TITV_LOW')?thr[r[3]]:(T[r[3]]!=null?T[r[3]]:thr[r[3]]);
      if(r[0]=='TITV_LOW'&&!(t>0))return; if(r[2]=='<'?v<t:v>t)f.push(r[0]);});
    if(!s.anc){var snp=s.m.snps; if(snp!=null&&SNPMED!=null){if(snp<SNPMED-thr.snp_z*SNPSIG)f.push('SNP_LOW');else if(snp>SNPMED+thr.snp_z*SNPSIG)f.push('SNP_HIGH');}
      if(s.m.het_frac!=null&&s.m.het_frac>thr.het_max_frac)f.push('HET_HIGH');
      // snpEff functional-annotation flags (WARN-only, cohort-relative robust-z; aDNA-guarded: deamination inflates missense/stop)
      if(ANNHI&&s.m.ann_high!=null&&s.m.ann_high>ANNHI.med+thr.snp_z*ANNHI.sig)f.push('HIGH_IMPACT_EXCESS');
      if(LOFHI&&s.m.lof_pct!=null&&s.m.lof_pct>LOFHI.med+thr.snp_z*LOFHI.sig)f.push('LOF_EXCESS');
      if(MSIL&&s.m.missense_silent!=null&&s.m.missense_silent>MSIL.med+thr.snp_z*MSIL.sig)f.push('PNPS_PROXY_HIGH');
      if(s.m.annotated_pct!=null&&(s.m.annotated_pct<ANNOT_FLOOR||(ANNPCT&&s.m.annotated_pct<90&&s.m.annotated_pct<ANNPCT.med-thr.snp_z*ANNPCT.sig)))f.push('ANNOTATION_POOR');}
    s.m.n_lineages=nlin(s); if(s.m.n_lineages!=null&&s.m.n_lineages>1)f.push('MIXED');
    if(s.anc&&s.dmg&&s.dmg.ct1!=null&&athr.damage_min_ct!=null&&s.dmg.ct1<athr.damage_min_ct)f.push('DAMAGE_LOW');
    pipelineFlags(s).forEach(function(x){if(f.indexOf(x)<0)f.push(x);});
    s.f=f; s.v=f.some(function(x){return FAILF[x];})?'FAIL':(f.length?'WARN':'PASS');});
  // multivariate QC-outlier flag (cohort-wide fit; WARN-only, never FAIL, never auto-excluded)
  applyMahal();
  if(R.mahal_cut!=null)R.samples.forEach(function(s){
    if(s.f.indexOf('NO_DATA')>=0)return;                                  // a produced-nothing FAIL is not a "combination outlier"
    if(s.m.qc_mahal!=null&&s.m.qc_mahal>R.mahal_cut){if(s.f.indexOf('QC_OUTLIER')<0)s.f.push('QC_OUTLIER');if(s.v=='PASS')s.v='WARN';}});
  R.counts={PASS:0,WARN:0,FAIL:0}; R.samples.forEach(function(s){R.counts[s.v]++;});
}
// flag margin ("why"): human-readable reason for one flag on a sample (for the modal + exclusion reason)
function flagWhy(s,fl){var T=actv(s),M=MET;function u(k){return (M[k]&&M[k].kind=='pct')?'%':'';}
  switch(fl){
    case 'LOW_DEPTH':return 'depth '+shortv(s.m.mean_depth,'float')+' < '+T.depth_min;
    case 'LOW_BREADTH':return 'breadth '+shortv(s.m.breadth_pct,'pct')+'% < '+T.breadth_min+'%';
    case 'HIGH_MISSING':return 'missing '+shortv(s.m.missing_pct,'pct')+'% > '+T.missing_max+'%';
    case 'MAPPING_LOW':return 'mapped '+shortv(s.m.mapped_pct,'pct')+'% < '+T.mapping_min+'%';
    case 'HIGH_DUP':return 'dup '+shortv(s.m.duplication_pct,'pct')+'% > '+T.dup_max+'%';
    case 'HIGH_IUPAC':return 'IUPAC '+shortv(s.m.iupac_pct,'pct')+'% > '+T.iupac_max+'%';
    case 'TITV_LOW':return 'Ti/Tv '+shortv(s.m.ti_tv,'float')+' < '+thr.titv_min;
    case 'SNP_LOW':return 'SNPs '+shortv(s.m.snps,'int')+' below cohort (z<-'+thr.snp_z+')';
    case 'SNP_HIGH':return 'SNPs '+shortv(s.m.snps,'int')+' above cohort (z>'+thr.snp_z+')';
    case 'HET_HIGH':return 'het fraction '+shortv(s.m.het_frac,'pct')+'% > '+thr.het_max_frac+'%';
    case 'MIXED':return s.m.n_lineages+' lineages above '+thr.mixed_min_frac+'%';
    case 'DAMAGE_LOW':return "5' C>T "+(s.dmg&&s.dmg.ct1!=null?(s.dmg.ct1*100).toFixed(1)+'%':'NA')+' < '+(athr.damage_min_ct*100).toFixed(0)+'% (aDNA auth)';
    case 'NO_DATA':return 'no consensus / mapping output produced';
    case 'QC_OUTLIER':return 'unusual combination of QC metrics - Mahalanobis d² '+shortv(s.m.qc_mahal,'float')+' > cohort robust cut '+shortv(R.mahal_cut,'float')+(s._mdrv?' (driven by '+s._mdrv.map(function(d){return d.label+' '+(d.z>=0?'+':'')+d.z.toFixed(1)+'σ';}).join(', ')+')':'')+'; WARN only, no single metric fails';
    case 'HIGH_IMPACT_EXCESS':return 'HIGH-impact variants '+shortv(s.m.ann_high,'int')+' above cohort (robust z > '+thr.snp_z+'); possible contamination, misassembly, or a wrong reference annotation, or a genuinely reduced/divergent genome, not confirmed';
    case 'LOF_EXCESS':return 'loss-of-function fraction '+shortv(s.m.lof_pct,'pct')+'% above cohort (robust z > '+thr.snp_z+'); possible frameshift storm from indel-calling / assembly / annotation, or real pseudogenisation, not confirmed';
    case 'PNPS_PROXY_HIGH':return 'missense/silent ratio '+shortv(s.m.missense_silent,'float')+' above cohort (robust z > '+thr.snp_z+'); a spectrum proxy only, possible base-call error or contamination; not dN/dS, see the pN/pS panel';
    case 'ANNOTATION_POOR':return 'annotated '+shortv(s.m.annotated_pct,'pct')+'% below floor '+ANNOT_FLOOR+' or below cohort; the GFF-to-snpEff database may be mismatched for this reference';
    case 'LINEAGE_MISMATCH':return 'types as '+linMain(s.lineage)+(s.ref_lin?(' while the other samples on '+(s.ref?'reference '+s.ref:'its reference')+' type as '+linMain(s.ref_lin)):' unlike the other samples on its reference')+': probably mapped to a genome it does not belong to, so its SNPs measure that distance';
    default:return flagLab(fl);}}
// The top level of a lineage call ('A4;M_bovis;A4.7;A4.7.1' -> 'A4'), which is what a mismatch compares.
function linMain(l){return l?String(l).split(';')[0]:'NA';}
// The metric behind each flag: the cell of the sample table that tripped it is the one coloured.
var FLAGMET={LOW_DEPTH:'mean_depth',LOW_BREADTH:'breadth_pct',HIGH_MISSING:'missing_pct',MAPPING_LOW:'mapped_pct',
  HIGH_DUP:'duplication_pct',HIGH_IUPAC:'iupac_pct',TITV_LOW:'ti_tv',SNP_LOW:'snps',SNP_HIGH:'snps',HET_HIGH:'het_frac',
  MIXED:'n_lineages',HIGH_IMPACT_EXCESS:'ann_high',LOF_EXCESS:'lof_pct',PNPS_PROXY_HIGH:'missense_silent',
  ANNOTATION_POOR:'annotated_pct',QC_OUTLIER:'qc_mahal'};
// A flag as a reader needs it: [what is wrong, the value against the rule it broke]. flagWhy() keeps
// the terse form the exported exclusion.tsv has always carried.
function flagReason(s,fl){var T=actv(s),m=s.m;function v(k){return nv(m[k]);}
  switch(fl){
    case 'LOW_DEPTH':return [flagLab(fl),v('mean_depth')+'&#215;, needs '+T.depth_min+'&#215;'];
    case 'LOW_BREADTH':return [flagLab(fl),v('breadth_pct')+'% of the genome covered, needs '+T.breadth_min+'%'];
    case 'HIGH_MISSING':return [flagLab(fl),v('missing_pct')+'% of the consensus missing, at most '+T.missing_max+'%'];
    case 'MAPPING_LOW':return [flagLab(fl),v('mapped_pct')+'% of reads map, needs '+T.mapping_min+'%'];
    case 'HIGH_DUP':return [flagLab(fl),v('duplication_pct')+'% duplicates, at most '+T.dup_max+'%'];
    case 'HIGH_IUPAC':return [flagLab(fl),v('iupac_pct')+'% ambiguous bases, at most '+T.iupac_max+'%'];
    case 'TITV_LOW':return [flagLab(fl),'Ti/Tv '+v('ti_tv')+', needs '+thr.titv_min];
    case 'SNP_HIGH':case 'SNP_LOW':return [flagLab(fl),shortv(m.snps,'int')+' SNPs; the cohort median is '+shortv(SNPMED,'int')];
    case 'HET_HIGH':return [flagLab(fl),v('het_frac')+'% of variant calls heterozygous, at most '+thr.het_max_frac+'%'];
    case 'MIXED':return [flagLab(fl),m.n_lineages+' lineages above '+thr.mixed_min_frac+'% of the markers'];
    case 'LINEAGE_MISMATCH':return [flagLab(fl),'types as '+esc(linMain(s.lineage))+(s.ref_lin?'; the rest of '+esc(s.ref||'its reference')+' types as '+esc(linMain(s.ref_lin)):'')];
    case 'QC_OUTLIER':return [flagLab(fl),(s._mdrv&&s._mdrv.length?s._mdrv.slice(0,3).map(function(d){return esc(d.label)+' '+(d.z>=0?'+':'&#8722;')+Math.abs(d.z).toFixed(1)+'&#963;';}).join(', '):'d&#178; '+v('qc_mahal'))+' against the cohort'];
    case 'HIGH_IMPACT_EXCESS':return [flagLab(fl),shortv(m.ann_high,'int')+' HIGH-impact variants; the cohort median is '+(ANNHI?shortv(ANNHI.med,'int'):'NA')];
    case 'LOF_EXCESS':return [flagLab(fl),v('lof_pct')+'% of coding changes; the cohort median is '+(LOFHI?shortv(LOFHI.med,'float'):'NA')+'%'];
    case 'PNPS_PROXY_HIGH':return [flagLab(fl),'missense/silent '+v('missense_silent')+'; the cohort median is '+(MSIL?shortv(MSIL.med,'float'):'NA')];
    case 'ANNOTATION_POOR':return [flagLab(fl),'only '+v('annotated_pct')+'% of variants annotated'];
    case 'DAMAGE_LOW':return [flagLab(fl),"5&#8242; C&gt;T "+(s.dmg&&s.dmg.ct1!=null?(s.dmg.ct1*100).toFixed(1)+'%':'NA')+', needs '+(athr.damage_min_ct*100).toFixed(0)+'%'];
    case 'NO_DATA':return [flagLab(fl),'no mapping or consensus output'];
    default:return [flagLab(fl),esc(flagWhy(s,fl))];}}
