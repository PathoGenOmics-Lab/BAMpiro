// ---- pages: the report reads as eight pages, and the sidebar is both the page list and the contents of
// the page on screen. Only that page is redrawn on a change (a new threshold, a highlighted sample, the
// theme); the others are marked stale and redraw when opened, which is also when they know their width.
var PAGES=['summary','qc','genome','related','variants','drug','gconv','diag'];
var PAGE_RENDER={
  summary:function(){renderExec();},
  qc:function(){renderOverview();renderFlags();renderCuration();renderTable();renderLineages();renderKraken();renderPlots();renderADNA();},
  genome:function(){renderStacks();renderGenome();renderDeletions();renderFunction();renderGeneBurden();renderHotspots();renderPnps();},
  related:function(){renderRelatedness();},
  variants:function(){renderDynamics();renderEpistasis();renderSnpMatrix();renderVarDose();},
  drug:function(){renderDrug();},
  gconv:function(){renderGconv();},
  diag:function(){renderScatter();renderCorr();renderQCspace();renderRefBias();renderDoseTx();renderTemporal();}};
var curPage='summary', pageDirty={};
function renderPage(p){if(PAGE_RENDER[p])PAGE_RENDER[p]();pageDirty[p]=false;}
function renderAll(){renderNavBadges();PAGES.forEach(function(p){pageDirty[p]=true;});renderPage(curPage);renderInsights();}
function renderEverything(){renderNavBadges();PAGES.forEach(renderPage);renderInsights();}
function pageOf(id){var e=el(id);if(!e)return null;var pg=e.closest?e.closest('.page'):null;return pg?pg.getAttribute('data-page'):null;}
function pageShown(p){var pg=el('p-'+p);return !!pg&&pg.getAttribute('data-empty')!=='1';}
function showPage(p,target){
  if(!pageShown(p))p='summary';
  var changed=(p!==curPage); curPage=p;
  PAGES.forEach(function(q){var pg=el('p-'+q);if(pg)pg.classList.toggle('on',q===p);
    var g=document.querySelector('#toc .toc-group[data-page="'+q+'"]');if(g)g.classList.toggle('open',q===p);});
  if(changed||pageDirty[p]){renderPage(p);renderInsights();}
  var t=(target&&target.indexOf('p-')!==0)?el(target):null;
  if(t&&t.scrollIntoView)t.scrollIntoView(); else window.scrollTo(0,0);
  if(window.__spy)window.__spy();
}
// '#p-qc' opens a page, '#flagged' the page holding that section, scrolled to it; false if neither.
function goHash(h){var id=(h||'').replace(/^#/,'');if(!id){showPage('summary');return true;}
  var p=(id.indexOf('p-')===0)?id.slice(2):pageOf(id); if(!p||!PAGE_RENDER[p])return false; showPage(p,id); return true;}

// ---- static wiring ----
el('meta').textContent=R.samples.length+' samples · '+R.generated;
(function(){var hv=el('hver'); if(hv){ if(R.version){hv.textContent='v'+R.version; hv.title='BAMpiro version '+R.version;} else hv.style.display='none'; }})();
fillIcons();   // swap every static data-ic placeholder (header, TOC chevrons, buttons, modal) for its inline SVG
// click an (i) info icon -> show a persistent popover with its definition (capture phase so it
// beats the column-sort handler); click anywhere / Esc / scroll to dismiss.
(function(){
  var pop=el('infopop');
  document.addEventListener('click',function(e){
    var ic=(e.target&&e.target.closest)?e.target.closest('.infoi'):null;
    if(ic){
      e.stopPropagation(); e.preventDefault();
      if(pop._for===ic&&pop.style.display==='block'){ pop.style.display='none'; pop._for=null; return; }
      var _h=ic.closest('[title]'); pop.textContent=ic.getAttribute('data-info')||(_h?_h.getAttribute('title'):'')||'';
      var pw=Math.min(320,window.innerWidth-24); pop.style.maxWidth=pw+'px'; pop.style.display='block';
      var r=ic.getBoundingClientRect();
      var left=Math.min(Math.max(8,r.left-4),window.innerWidth-pw-8), top=r.bottom+8;
      if(top+pop.offsetHeight+8>window.innerHeight){ var up=r.top-8-pop.offsetHeight; if(up>4)top=up; }
      pop.style.left=left+'px'; pop.style.top=top+'px'; pop._for=ic;
    } else if(pop.style.display==='block'){ pop.style.display='none'; pop._for=null; }
  },true);
  document.addEventListener('keydown',function(e){ if(e.key==='Escape'&&pop.style.display==='block'){pop.style.display='none';pop._for=null;} });
  window.addEventListener('scroll',function(){ if(pop.style.display==='block'){pop.style.display='none';pop._for=null;} },true);
})();
el('foot').innerHTML='<span class="foot-brand">BAMpiro'+(R.version?' <b>v'+esc(R.version)+'</b>':'')+' · <a href="'+esc(R.repo_url)+'" target="_blank" rel="noopener noreferrer">'+icon('github','sort')+'source on GitHub'+icon('ext','sort')+'</a></span> · Generated '+R.generated+' · the thresholds can be changed on the Sample QC page; the pipeline gate and qc_flags.tsv keep the defaults ('+
  Object.keys(R.thresholds).map(function(k){return k+'='+R.thresholds[k];}).join(', ')+'). NA = not reported.';
el('colmenu').innerHTML='<div style="display:flex;gap:12px;margin-bottom:6px;padding-bottom:6px;border-bottom:1px solid var(--line);font-size:12px"><a href="#" id="colall" style="color:var(--accent)">show all</a><a href="#" id="colnone" style="color:var(--accent)">hide all</a></div>'+R.metrics.map(function(m){return '<label><input type="checkbox" data-k="'+m.key+'"'+(st.hidden[m.key]?'':' checked')+'> '+esc(m.label)+'</label>';}).join('');
Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.onchange=function(){if(cb.checked)delete st.hidden[cb.getAttribute('data-k')];else st.hidden[cb.getAttribute('data-k')]=1;renderTable();saveState();};});
(function(){var ca=el('colall'),cn=el('colnone');
  if(ca)ca.onclick=function(e){e.preventDefault();st.hidden={};Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.checked=true;});renderTable();saveState();};
  if(cn)cn.onclick=function(e){e.preventDefault();Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.checked=false;st.hidden[cb.getAttribute('data-k')]=1;});renderTable();saveState();};})();
el('q').oninput=function(e){st.q=e.target.value.toLowerCase().trim();renderTable();clearTimeout(_qdb);_qdb=setTimeout(function(){renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderStacks();renderGenome();renderFunction();renderTemporal();renderGconv();},160);};
// per-panel gene search (Functional gene burden / Variable genes / pN-pS): filter each gene table by gene name
[['gbq','gbq',renderGeneBurden],['hotq','hotq',renderHotspots],['delq','delq',renderDeletions],['pnpsq','pnpsq',renderPnps],['vardoseq','vardoseq',renderVarDose]].forEach(function(w){var inp=el(w[0]);if(inp)inp.oninput=function(e){st[w[1]]=e.target.value.trim();w[2]();};});
el('of').onchange=function(e){st.onlyFlagged=e.target.checked;renderAll();};
Array.prototype.forEach.call(document.querySelectorAll('#ptype button'),function(b){b.onclick=function(){st.ptype=b.getAttribute('data-t');
  Array.prototype.forEach.call(document.querySelectorAll('#ptype button'),function(x){x.classList.toggle('on',x==b);});renderPlots();};});
// scatter axis pickers
['sx','sy'].forEach(function(ax){var sel=el(ax);sel.innerHTML=R.metrics.map(function(m){return '<option value="'+m.key+'"'+(st[ax]==m.key?' selected':'')+'>'+esc(m.label)+'</option>';}).join('');
  sel.onchange=function(){st[ax]=sel.value;renderScatter();};});
// live thresholds + presets
var THL=[['depth_min','Depth min'],['breadth_min','Breadth min %'],['missing_max','Missing max %'],['mapping_min','Mapped min %'],['dup_max','Dup max %'],['iupac_max','IUPAC max %'],['titv_min','Ti/Tv min'],['snp_z','SNP z'],['het_max_frac','Het % max'],['mixed_min_frac','Mixed lin % min']];
var PRESETS=[
  {id:'gate',label:'gate defaults',th:assign({},R.thresholds),note:'The cut-offs the Snakemake qc_gate uses (config report_* keys).'},
  {id:'strict',label:'strict (modern WGS)',th:{depth_min:20,breadth_min:95,missing_max:5,mapping_min:90,dup_max:30,iupac_max:2,titv_min:1.5,snp_z:3,het_max_frac:1.5,mixed_min_frac:2},note:'Confident modern Illumina isolate: 20x, 95% breadth, <5% missing, Ti/Tv >=1.5.'},
  {id:'lenient',label:'lenient (aDNA / low-cov)',th:{depth_min:3,breadth_min:60,missing_max:40,mapping_min:50,dup_max:80,iupac_max:5,titv_min:0,snp_z:4,het_max_frac:8,mixed_min_frac:5},note:'Degraded / low-coverage library at the 3x calling floor: rescues calibration tips.'}
];
function applyThr(next){Object.keys(next).forEach(function(k){if(k in thr)thr[k]=next[k];});
  Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in thr)inp.value=thr[k];});
  recompute();renderAll();saveState();}
el('thbox').innerHTML='<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">preset<select id="thpreset" class="msel" style="padding:4px 6px"><option value="">custom…</option>'+
  PRESETS.map(function(p){return '<option value="'+p.id+'">'+esc(p.label)+'</option>';}).join('')+'</select></label>'+
  '<span id="thnote" style="flex-basis:100%;font-size:10.5px;color:#8895a6;margin-top:-2px"></span>'+
  THL.map(function(t){return '<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">'+t[1]+
  '<input type="number" step="any" data-t="'+t[0]+'" value="'+thr[t[0]]+'" style="width:78px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;font-size:12px"></label>';}).join('')+
  '<button class="btn" id="threset" style="align-self:flex-end">reset</button>';
var thpre=el('thpreset'),thnote=el('thnote');
if(thpre)thpre.onchange=function(){var p=null;PRESETS.forEach(function(x){if(x.id==thpre.value)p=x;});if(!p){thnote.textContent='';return;}thnote.textContent=p.note;applyThr(p.th);};
Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){inp.oninput=function(){var v=parseFloat(inp.value);if(!isNaN(v)){thr[inp.getAttribute('data-t')]=v;if(thpre)thpre.value='';if(thnote)thnote.textContent='';clearTimeout(_thdb);_thdb=setTimeout(function(){recompute();renderAll();saveState();},180);}};});
el('threset').onclick=function(){if(thpre)thpre.value='';if(thnote)thnote.textContent='';applyThr(assign({},R.thresholds));};
// ancient (aDNA) live thresholds
if(R.n_ancient){var ATH=[['depth_min','aDNA depth min'],['breadth_min','aDNA breadth %'],['missing_max','aDNA missing %'],['mapping_min','aDNA mapped %'],['dup_max','aDNA dup %'],['iupac_max','aDNA IUPAC %'],['damage_min_ct',"5′ C>T min (0-1)"]];
  el('athbox').innerHTML='<div style="flex-basis:100%;font-size:10px;color:#8a5a12;font-weight:600;text-transform:uppercase;letter-spacing:.06em">Ancient (aDNA) thresholds &middot; a 5x mummy is judged here, not against the modern gate</div>'+
    ATH.map(function(t){return '<label style="display:inline-flex;flex-direction:column;font-size:10px;color:#64748b;gap:2px">'+t[1]+
    '<input type="number" step="any" data-t="'+t[0]+'" value="'+(athr[t[0]]!=null?athr[t[0]]:'')+'" style="width:78px;padding:4px 6px;border:1px solid var(--line);border-radius:6px;font-size:12px"></label>';}).join('');
  Array.prototype.forEach.call(document.querySelectorAll('#athbox input'),function(inp){inp.oninput=function(){var v=parseFloat(inp.value);if(!isNaN(v)){athr[inp.getAttribute('data-t')]=v;clearTimeout(_thdb);_thdb=setTimeout(function(){recompute();renderAll();saveState();},180);}};});}
// CSV export of the current (filtered rows, visible columns) table
el('csv').onclick=function(){var mets=tableMetrics();
  var head=['sample','verdict'].concat(mets.map(function(m){return m.key;})).concat(['lineage','flags']);
  var lines=[head.join('\t')]; visible().filter(colMatch).forEach(function(s){lines.push([s.s,s.v].concat(mets.map(function(m){return s.m[m.key]==null?'':s.m[m.key];})).concat([s.lineage||'',s.f.join(';')]).join('\t'));});
  var blob=new Blob([lines.join('\n')],{type:'text/tab-separated-values'}),a=document.createElement('a');
  a.href=URL.createObjectURL(blob);a.download='qc_table.tsv';a.click();URL.revokeObjectURL(a.href);toast('Saved qc_table.tsv ('+visible().filter(colMatch).length+' rows)');};
// tooltips + click on plots/scatter
function bandTip(pk,val){var b=bandFor(pk,thr);if(!b)return'';var mk=(MET[pk]||{}).kind,inb=val>=b[0]&&val<=b[1];
  var lab=(b[0]==-Infinity)?('≤ '+shortv(b[1],mk)):(b[1]==Infinity)?('≥ '+shortv(b[0],mk)):(shortv(b[0],mk)+' to '+shortv(b[1],mk));
  return '<br><span style="color:'+(inb?'#8fe3c0':'#f2b8bc')+'">band '+lab+(inb?' ✓ in':' ✗ out')+'</span>';}
function wireHover(host){host.addEventListener('mousemove',function(e){var t=e.target,lin;
  if(t.tagName=='circle'&&t.hasAttribute('data-val')){lin=t.getAttribute('data-lin');var pk=t.getAttribute('data-pk'),vv=+t.getAttribute('data-val');
    tip('<b>'+esc(t.getAttribute('data-s'))+'</b>'+(lin?' · '+esc(lin):'')+'<br>'+esc(t.getAttribute('data-lab'))+': '+shortv(vv,t.getAttribute('data-kind'))+(pk?bandTip(pk,vv):''),e.clientX,e.clientY);}
  else if(t.tagName=='rect'&&t.hasAttribute('data-val')){lin=t.getAttribute('data-lin');var pk2=t.getAttribute('data-pk'),vv2=+t.getAttribute('data-val');
    tip('<b>'+esc(t.getAttribute('data-s'))+'</b>'+(lin?' · '+esc(lin):'')+'<br>'+esc(t.getAttribute('data-lab'))+': '+shortv(vv2,t.getAttribute('data-kind'))+(pk2?bandTip(pk2,vv2):''),e.clientX,e.clientY);}
  else if(t.tagName=='circle'&&t.hasAttribute('data-x'))tip('<b>'+esc(t.getAttribute('data-s'))+'</b>'+(t.getAttribute('data-lin')?' · '+esc(t.getAttribute('data-lin')):'')+'<br>'+esc(t.getAttribute('data-xl'))+': '+shortv(+t.getAttribute('data-x'),t.getAttribute('data-xk'))+'<br>'+esc(t.getAttribute('data-yl'))+': '+shortv(+t.getAttribute('data-y'),t.getAttribute('data-yk')),e.clientX,e.clientY);
  else tip('');});
  host.addEventListener('mouseleave',function(){tip('');});
  host.addEventListener('click',function(e){if((e.target.tagName=='circle'||e.target.tagName=='rect')&&e.target.getAttribute('data-s'))setHi(e.target.getAttribute('data-s'));});}
wireHover(el('plots')); wireHover(el('scatter')); wireHover(el('qcpca_body')); wireHover(el('divcomp_body')); wireHover(el('temporal_body'));
// correlation matrix: hover tip + click a cell -> load that metric pair into the scatter
(function(){var c=el('corr_body'); if(!c)return;
  c.addEventListener('mousemove',function(e){var t=e.target;
    if(t.tagName=='rect'&&t.hasAttribute('data-r')&&t.getAttribute('data-xk')){var r=t.getAttribute('data-r');
      tip('<b>'+esc((MET[t.getAttribute('data-yk')]||{}).label||'')+'</b> vs <b>'+esc((MET[t.getAttribute('data-xk')]||{}).label||'')+'</b><br>Spearman rho '+(r===''?'NA':r),e.clientX,e.clientY);}
    else tip('');});
  c.addEventListener('mouseleave',function(){tip('');});
  c.addEventListener('click',function(e){var t=e.target;
    if(t.tagName=='rect'&&t.getAttribute('data-xk')&&t.getAttribute('data-xk')!=t.getAttribute('data-yk')){
      st.sx=t.getAttribute('data-xk');st.sy=t.getAttribute('data-yk');
      var sxs=el('sx'),sys=el('sy'); if(sxs)sxs.value=st.sx; if(sys)sys.value=st.sy;
      renderScatter();el('corr').scrollIntoView({behavior:'smooth'});}});})();
// scatter rubber-band select -> add the enclosed samples to the exclusion basket (rAF-gated overlay, hit-test on release)
(function(){var host=el('scatter'); if(!host)return;
  var dragging=false,rectL=0,rectT=0,x0=0,y0=0,x1=0,y1=0,raf=0;
  function paint(){raf=0;var ov=el('scbrush');if(!ov)return;
    ov.setAttribute('x',Math.min(x0,x1));ov.setAttribute('y',Math.min(y0,y1));ov.setAttribute('width',Math.abs(x1-x0));ov.setAttribute('height',Math.abs(y1-y0));ov.style.display='';}
  function move(e){x1=e.clientX-rectL;y1=e.clientY-rectT;if(!raf)raf=requestAnimationFrame(paint);}
  function up(){if(!dragging)return;dragging=false;
    document.removeEventListener('mousemove',move);document.removeEventListener('mouseup',up);
    var ov=el('scbrush');if(ov)ov.style.display='none';
    if(!SGEO||(Math.abs(x1-x0)<4&&Math.abs(y1-y0)<4))return;   // a click, not a drag
    var G=SGEO,l=Math.min(x0,x1),rr=Math.max(x0,x1),tp=Math.min(y0,y1),bt=Math.max(y0,y1);
    function sx(v){return G.pad+(G.xr[1]>G.xr[0]?(v-G.xr[0])/(G.xr[1]-G.xr[0]):0.5)*G.plot;}
    function sy(v){return G.H-G.pad-(G.yr[1]>G.yr[0]?(v-G.yr[0])/(G.yr[1]-G.yr[0]):0.5)*G.ph;}
    function inbox(s){var xv=s.m[G.xk],yv=s.m[G.yk];if(xv==null||yv==null)return false;
      var px=sx(xv),py=sy(yv);return px>=l&&px<=rr&&py>=tp&&py<=bt;}
    var vis=visible(),visSet={}; vis.forEach(function(s){visSet[s.s]=1;});
    var picked=vis.filter(inbox);
    var skipped=R.samples.filter(function(s){return !visSet[s.s]&&inbox(s);}).length;   // dimmed, out-of-filter dots inside the box
    var dim=skipped?' <span style="color:#94a3b8">('+skipped+' filtered-out dot'+(skipped>1?'s':'')+' ignored)</span>':'';
    var rd=el('scbrushinfo'); if(!picked.length){if(rd)rd.innerHTML=skipped?'no in-view samples in the box'+dim:'';return;}
    var added=picked.filter(function(s){return !st.excl[s.s];});   // undo only reverts what this brush added
    added.forEach(function(s){st.excl[s.s]=1;}); renderTable();renderCuration();
    if(rd){rd.innerHTML='basketed <b>'+added.length+'</b> of '+picked.length+' &#8594; exclusion'+dim+(added.length?' <button class="btn" id="scbrushundo" style="padding:2px 8px;font-size:11px">undo</button>':'');
      var ub=el('scbrushundo'); if(ub)ub.onclick=function(){added.forEach(function(s){delete st.excl[s.s];});renderTable();renderCuration();rd.innerHTML='';};}}
  host.addEventListener('mousedown',function(e){if(e.button!==0)return;var svg=el('scsvg');if(!svg)return;
    var rc=svg.getBoundingClientRect();rectL=rc.left;rectT=rc.top;
    dragging=true;x0=x1=e.clientX-rectL;y0=y1=e.clientY-rectT;
    document.addEventListener('mousemove',move);document.addEventListener('mouseup',up);e.preventDefault();});})();
// genome landscape hover + click
(function(){var g=el('genome_body'); if(!g)return;
  var LAB={missing:'missing',del:'deleted',snp:'homozygous SNPs',snpkb:'SNPs per callable kb',het:'het variants',indel:'indels'};
  g.addEventListener('mousemove',function(e){var t=e.target;
    if(t.tagName=='rect'&&t.hasAttribute('data-v')){var bin=+t.getAttribute('data-bin'),nb=R.nbins||200,gl=R.genome_len||nb,tk=st.gtrack||'missing';
      var p0=Math.round(bin/nb*gl),p1=Math.round((bin+1)/nb*gl),v=t.getAttribute('data-v');
      var val=(v==='')?'NA':(tk=='missing'?v+'% missing':tk=='del'?v+'% in a stretch without reads that other samples read':v+' '+LAB[tk]);
      tip('<b>'+esc(t.getAttribute('data-s')||'cohort (mean)')+'</b><br>'+fmtpos(p0)+' - '+fmtpos(p1)+'<br>'+val,e.clientX,e.clientY);}
    else tip('');});
  g.addEventListener('mouseleave',function(){tip('');});
  g.addEventListener('click',function(e){if(e.target.getAttribute('data-s'))setHi(e.target.getAttribute('data-s'));});})();
// genome brush: drag to ZOOM the plot into a reference span (and filter the Variable-genes table); nested drags zoom further
function genomeReadout(b0,b1,zoomed){var nb=R.nbins||200,gl=R.genome_len||nb;
  var p0=Math.round(b0/nb*gl),p1=Math.round((b1+1)/nb*gl);
  var ng=(R.genes||[]).filter(function(ge){return ge.end>=p0&&ge.start<=p1;}).length;
  var rd=el('gselreadout');if(rd)rd.innerHTML=(zoomed?'zoomed ':'')+'<b>'+fmtpos(p0)+' - '+fmtpos(p1)+'</b> &middot; '+(p1-p0).toLocaleString('en-US')+' bp'+(R.genes&&R.genes.length?' &middot; '+ng+' gene'+(ng==1?'':'s'):'')+' <button class="btn" id="gselclear" style="padding:2px 8px;font-size:11px">'+(zoomed?'reset zoom':'release to zoom')+'</button>';}
function genomeResetZoom(){st.gsel=null;st.gzoom=null;st.geneMark=null;var gg=el('genegoto');if(gg)gg.value='';renderGenome();renderHotspots();var rd=el('gselreadout');if(rd)rd.innerHTML='';}
(function(){var g=el('genome_body'); if(!g)return;
  var dragging=false,rectL=0,startBin=0,curBin=0,raf=0;
  function binAt(cx){if(!GGEO)return 0;var b=GGEO.z0+Math.floor((cx-rectL-GGEO.gut)/GGEO.plotW*GGEO.winN);return Math.max(GGEO.z0,Math.min(GGEO.z0+GGEO.winN-1,b));}
  function paint(){raf=0;var gb=el('gbrush');if(!gb||!GGEO)return;
    var b0=Math.min(startBin,curBin),b1=Math.max(startBin,curBin);
    var x0=GGEO.gut+(b0-GGEO.z0)/GGEO.winN*GGEO.plotW,x1=GGEO.gut+(b1+1-GGEO.z0)/GGEO.winN*GGEO.plotW;
    gb.setAttribute('x',x0.toFixed(1));gb.setAttribute('width',(x1-x0).toFixed(1));gb.style.display='';
    genomeReadout(b0,b1,false);}
  function move(e){curBin=binAt(e.clientX);if(!raf)raf=requestAnimationFrame(paint);}
  function up(){if(!dragging)return;dragging=false;
    document.removeEventListener('mousemove',move);document.removeEventListener('mouseup',up);
    var b0=Math.min(startBin,curBin),b1=Math.max(startBin,curBin);
    if(b1>b0){st.gsel={b0:b0,b1:b1};st.gzoom={b0:b0,b1:b1};renderGenome();renderHotspots();genomeReadout(b0,b1,true);
      var cb=el('gselclear');if(cb)cb.onclick=function(ev){ev.stopPropagation();genomeResetZoom();};}}
  g.addEventListener('mousedown',function(e){if(e.button!==0)return;var svg=g.querySelector('svg');if(!svg)return;
    rectL=svg.getBoundingClientRect().left;dragging=true;startBin=curBin=binAt(e.clientX);
    document.addEventListener('mousemove',move);document.addEventListener('mouseup',up);e.preventDefault();});})();
// go-to-gene: type a gene name -> zoom the plot to that gene and mark it (needs a GFF, i.e. R.genes)
(function(){var gg=el('genegoto');if(!gg)return;
  if(!(R.genes&&R.genes.length)){gg.style.display='none';return;}
  gg.oninput=function(e){var q=e.target.value.trim().toLowerCase();
    if(!q){genomeResetZoom();return;}
    var g=null,i;for(i=0;i<R.genes.length;i++){if((R.genes[i].name||'').toLowerCase().indexOf(q)>=0){g=R.genes[i];break;}}
    if(!g)return;
    var nb=R.nbins||200,gl=R.genome_len||nb;
    var b0=Math.max(0,Math.min(nb-1,Math.floor(g.start/gl*nb))),b1=Math.max(0,Math.min(nb-1,Math.floor(g.end/gl*nb)));
    var pad=Math.max(3,Math.round((b1-b0)*0.6)+2);
    st.gzoom={b0:Math.max(0,b0-pad),b1:Math.min(nb-1,b1+pad)}; st.geneMark={b0:b0,b1:b1,name:g.name};
    if(gtrackHas('snp')){st.gtrack='snp';Array.prototype.forEach.call(document.querySelectorAll('#gtrack button'),function(b){b.classList.toggle('on',b.getAttribute('data-gt')=='snp');});}
    renderGenome();renderHotspots();
    var rd=el('gselreadout');if(rd)rd.innerHTML='gene <b>'+esc(g.name)+'</b> &middot; '+fmtpos(g.start)+' - '+fmtpos(g.end)+' <button class="btn" id="gselclear" style="padding:2px 8px;font-size:11px">reset zoom</button>';
    var cb=el('gselclear');if(cb)cb.onclick=function(ev){ev.stopPropagation();genomeResetZoom();};};})();
if(!R.samples.some(function(s){return s.miss||s.trk;})){var gs=el('genome');if(gs)gs.style.display='none';var ng=el('nav-genome');if(ng)ng.style.display='none';}
if(!relSec()){['reldist','relclus','relgroup','nav-reldist','nav-relclus','nav-relgroup'].forEach(function(id){var x=el(id);if(x)x.style.display='none';});}
if(!(R.coverage&&R.coverage.regions)){var dlx=el('deletions');if(dlx)dlx.style.display='none';var ndlx=el('nav-del');if(ndlx)ndlx.style.display='none';}
if(!R.samples.some(function(s){return s.trk&&s.trk.snp;})){var hsx=el('hotspots');if(hsx)hsx.style.display='none';var nhx=el('nav-hot');if(nhx)nhx.style.display='none';}
// signature layer gates: hide a panel when its data type is absent cohort-wide (organism-agnostic; render fns also degrade per-view)
if(!pca2(R.samples).ok){var qp=el('qcpca');if(qp)qp.style.display='none';var nqp=el('nav-pca');if(nqp)nqp.style.display='none';}
if(!R.samples.some(function(s){return s.m.snps!=null||s.m.snp_density!=null;})){var dc=el('divcomp');if(dc)dc.style.display='none';var ndc=el('nav-divcomp');if(ndc)ndc.style.display='none';}
if(!R.samples.some(function(s){return yearOf(s.date)!=null;})){var tpx=el('temporal');if(tpx)tpx.style.display='none';var ntpx=el('nav-temporal');if(ntpx)ntpx.style.display='none';}
if(!R.samples.some(function(s){return s.m.ann_high!=null||s.m.ann_moderate!=null||s.m.ann_modifier!=null;})){var fnx=el('function');if(fnx)fnx.style.display='none';var nfx=el('nav-function');if(nfx)nfx.style.display='none';}
if(!(R.gene_burden&&R.gene_burden.length)){var gbx=el('geneburden');if(gbx)gbx.style.display='none';var ngb=el('nav-geneburden');if(ngb)ngb.style.display='none';}
if(!(R.pnps&&R.pnps.length)){var ppx=el('pnps');if(ppx)ppx.style.display='none';var npp=el('nav-pnps');if(npp)npp.style.display='none';}
// aDNA panel: only when there are ancient samples (no placeholder otherwise)
if(!R.n_ancient){var adx=el('adna');if(adx)adx.style.display='none';var nadx=el('nav-adna');if(nadx)nadx.style.display='none';}
// SNP dynamics: only when the metadata gave connected time-series (else no section at all)
if(!(R.dynamics&&R.dynamics.groups&&R.dynamics.groups.length)){var dyx=el('dynamics');if(dyx)dyx.style.display='none';var ndyx=el('nav-dyn');if(ndyx)ndyx.style.display='none';}
if(!(R.epistasis&&R.epistasis.pairs&&R.epistasis.pairs.length)){var epx=el('epistasis');if(epx)epx.style.display='none';var nepx=el('nav-epi');if(nepx)nepx.style.display='none';}
if(!(R.dr&&R.dr.calls&&R.dr.calls.length)){var drx=el('drug');if(drx)drx.style.display='none';var ndrx=el('nav-drug');if(ndrx)ndrx.style.display='none';}
if(!(R.gconv&&R.gconv.tracts&&R.gconv.tracts.length)){var gcx=el('gconv');if(gcx)gcx.style.display='none';var ngcx=el('nav-gconv');if(ngcx)ngcx.style.display='none';}
if(!(R.kraken&&R.kraken.samples&&R.kraken.samples.length)){var kkx=el('kraken');if(kkx)kkx.style.display='none';var nkkx=el('nav-kraken');if(nkkx)nkkx.style.display='none';}
// genome track selector (Missing / SNPs / Het / Indels) - only offer tracks that have data
(function(){var host=el('gtrack'); if(!host)return;var avail=GTRACKS.filter(function(g){return gtrackHas(g.k);});
  if(avail.length<=1){host.style.display='none';return;}
  host.innerHTML=avail.map(function(g){return '<button'+(g.k==st.gtrack?' class="on"':'')+' data-gt="'+g.k+'">'+g.lab+'</button>';}).join('');
  Array.prototype.forEach.call(host.querySelectorAll('button'),function(b){b.onclick=function(){st.gtrack=b.getAttribute('data-gt');
    Array.prototype.forEach.call(host.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});renderGenome();};});})();
// mask-regions toggle (mtbc_mask etc.): grey the masked zones + exclude them from the SNP density + variable-gene ranking
(function(){var mb=el('maskbtn'); if(!mb)return;
  if(!R.mask_bins){mb.style.display='none';return;}
  mb.title=(R.mask_pct||0)+'% of the reference masked (PE/PPE, IS, DR, repeats); toggle to exclude these zones';
  mb.onclick=function(){st.maskOn=!st.maskOn;mb.classList.toggle('on',st.maskOn);renderGenome();renderHotspots();};})();
var rz;window.addEventListener('resize',function(){clearTimeout(rz);rz=setTimeout(function(){renderPlots();renderScatter();renderCorr();renderQCspace();renderRefBias();renderGenome();renderTemporal();renderGconv();if(curPage=='related')renderRelatedness();},120);});
// metric help panel
el('helpmenu').innerHTML=R.metrics.map(function(m){var d=R.defs[m.key]||['',''];return '<div class="hitem"><b>'+esc(m.label)+'</b> <span class="hk">'+esc(m.key)+'</span><div class="hd">'+esc(d[0]||'')+(d[1]?' <span class="hr">('+esc(d[1])+')</span>':'')+'</div></div>';}).join('');
// ---- colour-by toggle (beeswarm + scatter); hides the whole lineage UI when there is no lineage data ----
(function(){var cb=el('colorby');
  if(!R.lin_present){if(cb)cb.style.display='none';['leg-lin','linsum','nav-lin','groupui'].forEach(function(id){var e=el(id);if(e)e.style.display='none';});return;}
  if(cb){Array.prototype.forEach.call(cb.querySelectorAll('button'),function(b){b.onclick=function(){st.colorBy=b.getAttribute('data-cb');
    Array.prototype.forEach.call(cb.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});
    el('leg-qc').style.display=(st.colorBy=='lineage')?'none':'flex';el('leg-lin').style.display=(st.colorBy=='lineage')?'flex':'none';
    var pc=el('pcacb');if(pc)Array.prototype.forEach.call(pc.querySelectorAll('button'),function(x){x.classList.toggle('on',x.getAttribute('data-cb')==st.colorBy);});
    renderPlots();renderScatter();renderQCspace();};});}
  el('leg-lin').innerHTML=R.lineages.map(function(l){return '<span><i style="background:'+linColor(l)+'"></i>'+esc(l)+'</span>';}).join('')+'<span style="margin-left:auto">beeswarm dashed line = median</span>';
})();
// PCA colour-by (mirrors #colorby into st.colorBy; works even with no lineage data by dropping only the lineage button)
(function(){var cb=el('pcacb');if(!cb)return;
  if(!R.lin_present){var lb=cb.querySelector('[data-cb="lineage"]');if(lb)lb.style.display='none';}
  Array.prototype.forEach.call(cb.querySelectorAll('button'),function(b){b.onclick=function(){st.colorBy=b.getAttribute('data-cb');
    Array.prototype.forEach.call(cb.querySelectorAll('button'),function(x){x.classList.toggle('on',x==b);});
    var main=el('colorby');if(main)Array.prototype.forEach.call(main.querySelectorAll('button'),function(x){x.classList.toggle('on',x.getAttribute('data-cb')==st.colorBy);});
    var lq=el('leg-qc'),ll=el('leg-lin');if(lq)lq.style.display=(st.colorBy=='lineage')?'none':'flex';if(ll)ll.style.display=(st.colorBy=='lineage')?'flex':'none';
    renderPlots();renderScatter();renderQCspace();};});})();
// group-by-lineage checkbox
var glin=el('glin'); if(glin)glin.onchange=function(){st.groupLin=glin.checked;renderTable();};
var colfCb=el('colf'); if(colfCb)colfCb.onchange=function(){st.showColF=colfCb.checked;renderTable();};
// ancient / modern filter
if(R.n_ancient){el('ancfilter').innerHTML='<span class="seg" id="ancseg"><button class="on" data-a="">all</button><button data-a="mod">modern</button><button data-a="anc">aDNA <span class="k">'+R.n_ancient+'</span></button></span>';
  Array.prototype.forEach.call(document.querySelectorAll('#ancseg button'),function(b){b.onclick=function(){st.ancOnly=b.getAttribute('data-a')||null;
    Array.prototype.forEach.call(document.querySelectorAll('#ancseg button'),function(x){x.classList.toggle('on',x==b);});renderAll();};});}
// samplesheet-metadata cohort filter: one dropdown per categorical annotation column (site,
// treatment, ...), restricting the whole report to an exact value. Time/group columns (the
// dynamics axes) and the lineage column (its own Lineages filter) are left out; every field
// still appears as a SNP-matrix header level regardless.
(function(){
  var host=el('metafilter'), meta=R.sample_meta; if(!host)return;
  if(!(meta&&meta.fields&&meta.fields.length))return;
  var skip={}; if(meta.time_field)skip[meta.time_field]=1; if(meta.group_field)skip[meta.group_field]=1;
  var has=function(o,k){return Object.prototype.hasOwnProperty.call(o,k);};   // avoid prototype-chain hits (toString, constructor, ...)
  // the lineage column is left out by NAME (its own Lineages filter owns it) and/or by value — the
  // samplesheet may label lineages differently from the QC-assigned ones (sub-lineages, "Beijing", ...).
  function isLinName(f){return /^(lineage|linaje|lineage_?id|sub_?lineage)$/i.test(f);}
  function isLinValues(f){var seen=false; for(var i=0;i<R.samples.length;i++){var v=(meta.rows[R.samples[i].s]||{})[f];
    if(v&&v!=='NA'&&v!=='.'&&v!=='-'){ if(!has(LINCOL,v))return false; seen=true; }} return seen;}
  var vals={}, flds=meta.fields.filter(function(f){
    if(skip[f]||isLinName(f)||isLinValues(f))return false;
    var seen={},list=[]; R.samples.forEach(function(s){var v=(meta.rows[s.s]||{})[f];
      if(v&&v!=='NA'&&v!=='.'&&v!=='-'&&!has(seen,v)){seen[v]=1;list.push(v);}});
    vals[f]=list.sort(); return list.length>=2&&list.length<=12;
  });
  if(!flds.length)return;
  host.innerHTML='<span class="metaflt" title="filter the whole report by a samplesheet annotation column">'+
    flds.map(function(f){return '<label class="metasel">'+esc(f)+' <select data-mf="'+esc(f)+'"><option value="">all</option>'+
      vals[f].map(function(v){return '<option value="'+esc(v)+'">'+esc(v)+'</option>';}).join('')+'</select></label>';}).join('')+'</span>';
  Array.prototype.forEach.call(host.querySelectorAll('select'),function(sel){sel.onchange=function(){ applyMetaFilter(sel.getAttribute('data-mf'), sel.value); };});
})();
// set/clear a cohort metadata filter and keep the toolbar dropdown in sync; shared by the dropdowns
// and the clickable Dose x treatment group labels.
function applyMetaFilter(field,val){
  if(val)st.metaFilter[field]=val; else delete st.metaFilter[field];
  Array.prototype.forEach.call(document.querySelectorAll('#metafilter select'),function(s){
    if(s.getAttribute('data-mf')===field){ s.value=val||''; s.classList.toggle('on',!!val); }});
  renderAll();
}
// provenance / run-manifest header (self-documenting for a citable exclusion set)
(function(){var p=R.provenance||{},items=[];
  items.push('reference '+(p.reference||'NA')+(R.genome_len?' ('+R.genome_len.toLocaleString('en-US')+' bp)':''));
  if(p.container&&p.container!='none')items.push('container '+p.container);
  if(p.commit)items.push('commit '+p.commit);
  items.push(R.samples.length+' samples'+(R.n_ancient?' · '+R.n_ancient+' aDNA':''));
  var pe=el('prov'); if(pe)pe.innerHTML=items.map(function(t){return '<span>'+esc(t)+'</span>';}).join('');})();
// print / save as PDF. The SNP matrix in "show all" mode only keeps the visible window in the DOM (with tall
// spacer rows), which would print as a few rows over a big blank; collapse it to the top sites for a clean
// printout, then restore. Run synchronously in the click handler AND via beforeprint (Ctrl+P) to be safe.
// Printing lays every page out one after another, drawn at the width they have on paper.
function _printPrep(){ if(document.body.classList.contains('print-all'))return;
  document.body.classList.add('print-all'); renderEverything();
  if(snpmxAll){ window.__printWasAll=true; snpmxAll=false; if(window.__snpmxDraw)window.__snpmxDraw(); if(window.__syncSitesBtn)window.__syncSitesBtn(); } }
function _printRestore(){ if(!document.body.classList.contains('print-all'))return;
  document.body.classList.remove('print-all');
  if(window.__printWasAll){ window.__printWasAll=false; snpmxAll=true; if(window.__snpmxDraw)window.__snpmxDraw(); if(window.__syncSitesBtn)window.__syncSitesBtn(); }
  renderAll(); }
var pbtn=el('printBtn'); if(pbtn)pbtn.onclick=function(){ _printPrep(); window.print(); _printRestore(); };
window.addEventListener('beforeprint',_printPrep);
window.addEventListener('afterprint',_printRestore);
// ---- persistence of the curated view (basket + thresholds), namespaced per sample-set so two reports don't bleed ----
var SKEY='bampiro_qc_v2:'+R.samples.length+':'+(R.samples[0]?R.samples[0].s:'')+':'+(R.samples.length?R.samples[R.samples.length-1].s:'');
function saveState(){try{localStorage.setItem(SKEY,JSON.stringify({v:3,excl:st.excl,thr:thr,athr:athr,hidden:st.hidden,sortKey:st.sortKey,asc:st.asc}));}catch(e){}}
function loadState(){try{var s=JSON.parse(localStorage.getItem(SKEY)||'null');if(!s)return false;
  if(s.thr)Object.keys(s.thr).forEach(function(k){if(k in thr)thr[k]=s.thr[k];});
  if(s.athr)Object.keys(s.athr).forEach(function(k){if(k in athr)athr[k]=s.athr[k];});
  if(s.excl&&typeof s.excl=='object'){var have={};R.samples.forEach(function(x){have[x.s]=1;});
    st.excl={};Object.keys(s.excl).forEach(function(k){if(have[k])st.excl[k]=1;});}  // drop unknown sample ids
  if(s.v>=3&&s.hidden&&typeof s.hidden=='object')st.hidden=s.hidden;   // restore the chosen columns (not an older report's, saved before the table opened on the key metrics)
  if(s.sortKey){st.sortKey=s.sortKey;st.asc=!!s.asc;}              // and the sort order (SKEY is per-cohort)
  return true;}catch(e){return false;}}
// ---- expand-to-fill (fullscreen within the window) for the big panels ----
function collapseExpanded(){var ex=document.querySelector('.panel.expanded');if(!ex)return;ex.classList.remove('expanded');document.body.classList.remove('has-expanded');
  Array.prototype.forEach.call(document.querySelectorAll('.exp-h'),function(b){b.innerHTML=icon('maximize')+'full';});
  renderGenome();renderPlots();renderScatter();renderTable();renderQCspace();renderRefBias();renderFunction();renderGeneBurden();renderHotspots();renderPnps();renderGconv();}
Array.prototype.forEach.call(document.querySelectorAll('.exp-h'),function(b){b.onclick=function(){
  var panel=el(b.getAttribute('data-panel')); if(!panel)return;
  var willExpand=!panel.classList.contains('expanded'); collapseExpanded();
  if(willExpand){panel.classList.add('expanded');document.body.classList.add('has-expanded');b.innerHTML=icon('minimize')+'close';}
  var rn=b.getAttribute('data-render');
  setTimeout(function(){if(rn=='genome')renderGenome();else if(rn=='plots')renderPlots();else if(rn=='scatter')renderScatter();else if(rn=='corr')renderCorr();else if(rn=='pca')renderQCspace();else if(rn=='divcomp')renderRefBias();else if(rn=='function')renderFunction();else if(rn=='geneburden')renderGeneBurden();else if(rn=='hotspots')renderHotspots();else if(rn=='deletions')renderDeletions();else if(rn=='pnps')renderPnps();else if(rn=='dosetx')renderDoseTx();else if(rn=='vardose')renderVarDose();else if(rn=='gconv')renderGconv();else if(rn=='table')renderTable();},20);};});
if(el('expClose'))el('expClose').onclick=collapseExpanded;
document.addEventListener('keydown',function(e){if(e.key=='Escape')collapseExpanded();});
// per-sample detail modal close (button, backdrop, Esc)
el('modalx').onclick=closeDetail; el('modal').addEventListener('click',function(e){if(e.target==el('modal'))closeDetail();});
el('modal').addEventListener('keydown',function(e){   // trap Tab focus inside the open dialog
  if(e.key!=='Tab')return;
  var f=Array.prototype.filter.call(el('modal').querySelectorAll('a[href],button:not([disabled]),input:not([disabled]),select,textarea,[tabindex]:not([tabindex="-1"])'),function(x){return x.offsetParent!==null;});
  if(!f.length)return; var first=f[0],last=f[f.length-1];
  if(e.shiftKey&&document.activeElement===first){e.preventDefault();last.focus();}
  else if(!e.shiftKey&&document.activeElement===last){e.preventDefault();first.focus();}
});
document.addEventListener('keydown',function(e){if(e.key=='Escape'&&st.detail)closeDetail();});

var hadSaved=loadState();
Array.prototype.forEach.call(document.querySelectorAll('#thbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in thr)inp.value=thr[k];});
if(R.n_ancient)Array.prototype.forEach.call(document.querySelectorAll('#athbox input'),function(inp){var k=inp.getAttribute('data-t');if(k in athr)inp.value=athr[k];});
Array.prototype.forEach.call(document.querySelectorAll('#colmenu input'),function(cb){cb.checked=!st.hidden[cb.getAttribute('data-k')];});  // sync the column checkboxes to any restored/hidden set
recompute();
if(!(hadSaved&&Object.keys(st.excl).length))R.samples.forEach(function(s){if(s.v=='FAIL')st.excl[s.s]=1;});  // preselect FAILs unless a saved basket exists
renderEverything();   // every page once, so each panel settles whether it has data at all
(function(){   // a page none of whose sections has data is dropped, with its sidebar entry; the rest are renumbered
  var n=0;
  PAGES.forEach(function(p){var pg=el('p-'+p); if(!pg)return;
    var any=(p==='summary')||Array.prototype.some.call(pg.querySelectorAll('section[id]'),function(x){return x.style.display!=='none';});
    pg.setAttribute('data-empty',any?'0':'1');
    var g=document.querySelector('#toc .toc-group[data-page="'+p+'"]'); if(g)g.style.display=any?'':'none';
    if(!any)return; n++;
    var num=g?g.querySelector('.toc-n'):null; if(num)num.textContent=n;
    var eb=pg.querySelector('.page-h .eyebrow'); if(eb)eb.textContent=n+' \u00b7 '+eb.textContent.replace(/^\d+\s*\u00b7\s*/,'');
    pageDirty[p]=true;});   // drawn off-screen at a guessed width: redraw on first open
})();
curPage=null; if(!goHash(location.hash))showPage('summary');
(function(){  // dark / light theme toggle (the early head script set the initial class from the saved pref or OS)
  var tb=el('themeToggle'); if(!tb)return;
  function setIcon(){tb.innerHTML=icon(isDark()?'sun':'moon');}
  setIcon();
  tb.onclick=function(){var d=!isDark();document.documentElement.classList.toggle('dark',d);
    try{localStorage.setItem('bampiro_theme',d?'dark':'light');}catch(e){}
    TH=d?TH_DARK:TH_LIGHT; setIcon(); renderAll();};
})();
(function(){  // analytical read-outs: one show/hide toggle for every .insight block across the report
  var b=el('insightToggle'); if(!b)return;
  var saved; try{ saved=localStorage.getItem('bampiro_insights'); }catch(e){}
  var on=(saved==null)?true:(saved!=='off');
  function apply(){ document.documentElement.classList.toggle('no-insight',!on); b.classList.toggle('on',on); b.setAttribute('aria-pressed',on?'true':'false'); }
  apply();
  b.onclick=function(){ on=!on; try{localStorage.setItem('bampiro_insights',on?'on':'off');}catch(e){} apply(); };
})();
(function(){  // header 'show all' shortcut -> expand every SNP view at once (trajectories, epistasis pairs, matrix sites) and jump to the dynamics
  var b=el('allSitesBtn'); if(!b)return;
  var hasAny=(R.snp_matrix&&R.snp_matrix.rows&&R.snp_matrix.rows.length)||(R.dynamics&&R.dynamics.groups&&R.dynamics.groups.length)||(R.epistasis&&R.epistasis.pairs&&R.epistasis.pairs.length);
  if(!hasAny){ b.style.display='none'; return; }
  b.style.display='';
  function anyOn(){ return dynState.showAll || epiState.showAll || (R.snp_matrix&&R.snp_matrix.rows&&snpmxAll); }
  window.__syncSitesBtn=function(){   // keep the master switch in step with the per-panel buttons
    var on=!!anyOn(); b.classList.toggle('on',on); var l=el('allSitesLbl'); if(l)l.textContent=on?'show less':'show all'; b.setAttribute('aria-pressed',on?'true':'false');
  };
  window.__syncSitesBtn();
  b.onclick=function(){   // master switch: if anything is expanded collapse everything, else expand everything
    var on=!anyOn();
    if(window.__dynSetAll) window.__dynSetAll(on);                                    // trajectories
    if(window.__epiSetAll) window.__epiSetAll(on);                                    // epistasis pairs
    if(R.snp_matrix&&R.snp_matrix.rows){ snpmxAll=on; var w=el('snpmxwrap'); if(w)w.scrollTop=0; if(window.__snpmxDraw)window.__snpmxDraw(); }   // matrix sites
    window.__syncSitesBtn();
    // deliberately no scroll: the header switch just flips the state, it doesn't navigate anywhere
  };
})();
(function(){  // left sidebar: collapse toggle (a drawer on a phone) and the scroll-spy of the page on screen
  var toc=el('toc'), tg=el('toc-toggle'); if(!toc||!tg)return;
  tg.onclick=function(){ document.body.classList.toggle('toc-collapsed'); };
  var scrim=el('tocscrim'); if(scrim)scrim.onclick=function(){ document.body.classList.add('toc-collapsed'); };  // tap outside the drawer to dismiss (mobile)
  if(window.innerWidth&&window.innerWidth<860) document.body.classList.add('toc-collapsed');   // start collapsed on small screens
  var links=Array.prototype.slice.call(toc.querySelectorAll('.toc-link'));
  function spy(){
    var pg=el('p-'+curPage), y=92, cur=null, best=-Infinity;   // the section whose top last passed under the header
    links.forEach(function(a){ a.classList.remove('active');
      var s=el(a.getAttribute('href').slice(1)); if(!s||!pg||!pg.contains(s)||s.style.display==='none'||a.style.display==='none')return;
      var top=s.getBoundingClientRect().top; if(top<=y&&top>best){best=top;cur=a;} });
    if(!cur){ for(var i=0;i<links.length;i++){ var s0=el(links[i].getAttribute('href').slice(1)); if(s0&&pg&&pg.contains(s0)&&s0.style.display!=='none'&&links[i].style.display!=='none'){cur=links[i];break;} } }
    if(cur) cur.classList.add('active');
  }
  window.__spy=spy;
  window.addEventListener('scroll',spy);
  window.addEventListener('resize',spy);
  spy();
})();
// every in-page link (sidebar, summary cards, "see the evidence") goes through the page router
document.addEventListener('click',function(e){
  var a=(e.target&&e.target.closest)?e.target.closest('a[href^="#"]'):null; if(!a)return;
  var h=a.getAttribute('href'); if(!h||h==='#')return;
  if(goHash(h)){ e.preventDefault(); try{history.replaceState(null,'',h);}catch(err){}
    if(window.innerWidth&&window.innerWidth<860)document.body.classList.add('toc-collapsed'); }
});
window.addEventListener('hashchange',function(){goHash(location.hash);});
(function(){  // drop a clickable (i) into every section heading; the explanation text comes from R.section_info
  var info=R.section_info||{};
  Array.prototype.forEach.call(document.querySelectorAll('section[id]'),function(sec){
    var txt=info[sec.id]; if(!txt) return;
    var h2=sec.querySelector('h2'); if(!h2||h2.querySelector('.sec-infoi')) return;
    var ic=document.createElement('span');
    ic.className='infoi sec-infoi'; ic.textContent='i';
    ic.setAttribute('data-info',txt);
    ic.setAttribute('role','button'); ic.setAttribute('aria-label','About this analysis');
    var cap=h2.querySelector('.c');
    if(cap) h2.insertBefore(ic,cap); else h2.appendChild(ic);
  });
})();
})();
