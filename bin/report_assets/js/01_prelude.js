
(function(){
var R=REPORT;
function assign(t,s){for(var _k in s){if(Object.prototype.hasOwnProperty.call(s,_k))t[_k]=s[_k];}return t;}  // ES5 shallow copy
function zeros(n){var _a=new Array(n);for(var _i=0;_i<n;_i++)_a[_i]=0;return _a;}                              // ES5 zero-filled array
var DIST=R.dist;
var VCOL={PASS:'#94a3b8',WARN:'#d97706',FAIL:'#dc2626'}, VFILL={PASS:'#16a34a',WARN:'#d97706',FAIL:'#dc2626'};
// Theme colours for the JS-drawn SVG panels (plots, heatmaps): swapped when the dark theme toggles.
function isDark(){return document.documentElement.classList.contains('dark');}
var TH_LIGHT={ink:'#15202e',mut:'#6a7889',panel:'#ffffff',soft:'#f6f8fb',line:'#e6ebf2',grid:'#e6ebf1',axis:'#8895a6',cellnull:'#e9edf2',track:'#edf1f6',hl:'#fff8e1',faint:'#5b6b7e'};
var TH_DARK ={ink:'#e7edf4',mut:'#94a4b6',panel:'#18232f',soft:'#1e2a38',line:'#2a3543',grid:'#2a3543',axis:'#6f7f92',cellnull:'#222e3c',track:'#232f3d',hl:'#33371c',faint:'#9fb0c2'};
var TH=isDark()?TH_DARK:TH_LIGHT;
var BAR={hi_good:'#22a06b',hi_bad:'#e0544f',neu:'#4f83c2'};
// ---- lineage palette: deterministic, colour-blind-safe, self-contained (Okabe-Ito + Tol-muted; golden-angle overflow)
var LINPAL_BASE=['#4477aa','#ee6677','#228833','#ccbb44','#66ccee','#aa3377','#e69f00','#0072b2','#d55e00','#009e73','#cc79a7','#882255'];
R.lin_colors=R.lin_colors||{};   // canonical palette (e.g. mycolorsTB) keyed by lineage label
var LINCOL={};(function(){(R.lineages||[]).forEach(function(lab,i){
  if(R.lin_colors[lab]){LINCOL[lab]=R.lin_colors[lab];}                       // canonical colour if provided
  else if(i<LINPAL_BASE.length){LINCOL[lab]=LINPAL_BASE[i];}                  // else colour-blind-safe fallback
  else{var h=(i*137.508)%360;LINCOL[lab]='hsl('+h.toFixed(0)+',58%,52%)';}});})();
function linColor(lab){return (lab!=null&&LINCOL[lab])?LINCOL[lab]:'#b8c2cf';}
R.gene_map=R.gene_map||{};
function geneRv(g){ return (g&&R.gene_map[g])||''; }   // Mycobrowser (H37Rv) locus tag for a gene, or ''
function geneRvTag(g){ var rv=geneRv(g); return rv?(' <a class="rvtag" href="https://mycobrowser.epfl.ch/genes/'+esc(rv)+'" target="_blank" rel="noopener" title="Mycobrowser locus tag of '+esc(g)+' (opens mycobrowser.epfl.ch)">'+esc(rv)+'</a>'):''; }
// amino-acid change in the used-reference numbering, plus the H37Rv/Mycobrowser one in brackets when it differs
var AA2LBL=R.aa2_label||'H37Rv';   // label for the canonical-reference amino-acid numbering
function aaDual(aa,aaH){ if(!aa) return ''; return esc(aa)+((aaH&&aaH!==aa)?(' <span class="aah37" title="same variant in the '+esc(AA2LBL)+' reference numbering">['+esc(AA2LBL)+' '+esc(aaH)+']</span>'):''); }
var REFNAME=R.ref_name||'';   // the reference the samples were mapped against, shown in the SNP-table headers
// mapping-reference position, plus the reference-of-interest (H37Rv) coordinate in brackets when the pipeline
// provides one that differs (no shared-coordinate assumption; pos2 is a 'contig:pos' string from the lifted VCF)
function refPos(pos,pos2){ var m=esc(''+pos); var p2=pos2?(''+pos2).split(':').pop():''; if(p2&&p2!==(''+pos)){ m+=' <span class="aah37" title="same variant in the '+esc(AA2LBL)+' reference coordinates">['+esc(AA2LBL)+' '+esc(p2)+']</span>'; } return m; }
var thr=assign({},R.thresholds);
var athr=assign({},R.anc_thresholds||{});      // ancient (aDNA) threshold view
function actv(s){return (s.anc&&R.n_ancient)?athr:thr;}   // active threshold set for a sample
R.defs=R.defs||{}; R.extra=R.extra||[]; R.provenance=R.provenance||{};
// ---- derived metrics registered as first-class metrics ----
function nlin(s){if(!s.linc)return null;var n=0,tot=0,k;for(k in s.linc)tot+=s.linc[k];if(tot<=0)return null;for(k in s.linc){if((s.linc[k]/tot)*100>=thr.mixed_min_frac&&s.linc[k]>=3)n++;}return n;}
R.samples.forEach(function(s){s.m.het_frac=(s.hetf!=null)?s.hetf*100:null;s.m.n_lineages=nlin(s);});
R.metrics.push({key:'het_frac',label:'Het %',kind:'pct',dir:'hi_bad'});
R.metrics.push({key:'n_lineages',label:'# lin',kind:'int',dir:'hi_bad'});
if(DIST.indexOf('het_frac')<0)DIST.push('het_frac');
// ---- SNP density per callable Mb (guarded: only for a plausible genome length) ----
if(R.snp_density_ok){var GMB=R.genome_len/1e6;R.samples.forEach(function(s){var snp=s.m.snps,cal=s.m.callable_pct;s.m.snp_density=(snp!=null&&cal!=null&&cal>0)?(snp/(cal/100*GMB)):null;});R.metrics.push({key:'snp_density',label:'SNP/Mb',kind:'float',dir:'neu'});if(DIST.indexOf('snp_density')<0)DIST.push('snp_density');}
// ===================== SIGNATURE ANALYSIS LAYER (QC-space PCA + robust Mahalanobis + divergence/completeness + temporal) =====================
R.defs=R.defs||{};
