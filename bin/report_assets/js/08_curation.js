function renderCuration(){var ex=nExcl(),keep=R.samples.length-ex;
  var nb=el('nbasket'); if(nb)nb.innerHTML=icon('basket')+ex+' in the exclusion list';   // always-visible toolbar mirror of the basket
  el('curation').innerHTML='<div class="cur-intro"><b>Exclusion list</b>: the samples you are dropping from the downstream analysis. Failing samples start in it; tick or untick any sample above or in the table below, then export the dropped samples with their reasons (<b>exclusion.tsv</b>) or the ones you keep (<b>keep_list.txt</b>).</div>'+
    '<div class="cur-read"><b>'+ex+'</b> to exclude <span class="arw">→</span> <b>'+keep+'</b> kept for downstream</div>'+
   '<div class="cur-btns"><button class="btn" data-cur="fail">exclude FAILs</button><button class="btn" data-cur="flagged">exclude all flagged</button>'+
   '<button class="btn" data-cur="clear">clear</button><button class="btn" data-cur="invert">invert (shown)</button>'+
   '<button class="btn prim" data-cur="excl">'+icon('download')+'exclusion.tsv</button><button class="btn prim" data-cur="keep">'+icon('download')+'keep_list.txt</button></div>';
  Array.prototype.forEach.call(el('curation').querySelectorAll('[data-cur]'),function(b){b.onclick=function(){curAction(b.getAttribute('data-cur'));};});
  // every change to the list passes through here, so the flagged list's boxes follow it too
  Array.prototype.forEach.call(document.querySelectorAll('#flagtable .fcb'),function(cb){cb.checked=!!st.excl[cb.getAttribute('data-s')];});
  if(typeof saveState=='function')saveState();}
function curAction(a){
  if(a=='fail')R.samples.forEach(function(s){if(s.v=='FAIL')st.excl[s.s]=1;});
  else if(a=='flagged')R.samples.forEach(function(s){if(s.v!='PASS')st.excl[s.s]=1;});
  else if(a=='clear')st.excl={};
  else if(a=='invert')visible().forEach(function(s){if(st.excl[s.s])delete st.excl[s.s];else st.excl[s.s]=1;});
  else if(a=='excl'){exportExcl();return;} else if(a=='keep'){exportKeep();return;}
  renderTable();renderCuration();}
function exportExcl(){var lines=['sample\tverdict\tancient\tflags\treason'];
  R.samples.filter(function(s){return st.excl[s.s];}).sort(function(a,b){return a.s.localeCompare(b.s);}).forEach(function(s){
    var reason=s.f.length?s.f.map(function(f){return flagWhy(s,f);}).join('; '):'manual_qc_exclusion';
    lines.push([s.s,s.v,(s.anc?'yes':'no'),(s.f.join(';')||'.'),reason].join('\t'));});
  dl(lines.join('\n')+'\n','exclusion.tsv','text/tab-separated-values');}
function exportKeep(){var lines=R.samples.filter(function(s){return !st.excl[s.s];}).map(function(s){return s.s;}).sort();
  dl(lines.join('\n')+'\n','keep_list.txt','text/plain');}
function dl(txt,name,type){var blob=new Blob([txt],{type:type}),a=document.createElement('a');a.href=URL.createObjectURL(blob);a.download=name;a.click();URL.revokeObjectURL(a.href);toast('Saved '+name);}
var _toastT;
function toast(msg){var e0=el('toast'); if(!e0)return; e0.textContent=msg; e0.className='show'; clearTimeout(_toastT); _toastT=setTimeout(function(){e0.className='';},2200);}
// keep a screen reader / keyboard user inside the open dialog: the page behind is made inert
function setBgInert(on){['toc','toc-toggle'].forEach(function(id){var e=el(id); if(e){if(on){e.setAttribute('inert','');e.setAttribute('aria-hidden','true');}else{e.removeAttribute('inert');e.removeAttribute('aria-hidden');}}});
  Array.prototype.forEach.call(document.querySelectorAll('header,.wrap'),function(e){if(on){e.setAttribute('inert','');e.setAttribute('aria-hidden','true');}else{e.removeAttribute('inert');e.removeAttribute('aria-hidden');}});}

// ---- per-sample detail modal ----
function pctRank(key,val){if(val==null)return null;var vs=R.samples.map(function(s){return s.m[key];}).filter(function(v){return v!=null;});if(!vs.length)return null;var b=0;vs.forEach(function(v){if(v<val)b++;});return Math.round(100*b/vs.length);}
function openDetail(sid){var s=null;R.samples.forEach(function(x){if(x.s==sid)s=x;});if(!s)return;st.detail=sid;
  var rows=R.metrics.map(function(m){var v=s.m[m.key],r=RANGES[m.key],nn=(v==null||r[1]<=r[0])?0:Math.max(0,Math.min(1,(v-r[0])/(r[1]-r[0]))),pr=pctRank(m.key,v);
    var _dd=(R.defs[m.key]||[''])[0];return '<div class="drow"><div class="dk" title="'+esc(_dd)+'">'+esc(m.label)+(_dd?'<span class="infoi" title="'+esc(_dd)+'">i</span>':'')+'</div>'+
      '<div class="dbarwrap"><div class="dbar" style="width:'+(nn*100).toFixed(1)+'%;background:'+BAR[m.dir]+'"></div></div>'+
      '<div class="dv">'+(v==null?'<span class="na">NA</span>':fmt(v,m.kind))+'</div><div class="dp">'+(pr==null?'':('p'+pr))+'</div></div>';}).join('');
  var lin='';if(s.linf){var ks=Object.keys(s.linf).sort(function(a,b){return s.linf[b]-s.linf[a];});
    lin='<div class="dsub">Lineage composition</div><div class="lcomp">'+ks.map(function(k){return '<div class="lrow"><span class="lk">'+esc(k)+'</span><div class="lbarw"><div class="lbar" style="width:'+(s.linf[k]*100).toFixed(1)+'%"></div></div><span class="lv">'+(s.linf[k]*100).toFixed(1)+'%</span></div>';}).join('')+'</div>';}
  var dmg='';if(s.anc){var d=s.dmg;
    if(d){dmg='<div class="dsub">aDNA damage (mapDamage2) - consistent with ancient DNA, not proof of authenticity</div>'+
      '<div class="drow"><div class="dk">5&#39; C&gt;T (pos 1)</div><div class="dbarwrap"><div class="dbar" style="width:'+Math.min(100,(d.ct1||0)*100/0.3).toFixed(1)+'%;background:'+(s.f.indexOf('DAMAGE_LOW')>=0?'#e0544f':'#22a06b')+'"></div></div><div class="dv">'+(d.ct1==null?'NA':(d.ct1*100).toFixed(1)+'%')+'</div><div class="dp"></div></div>'+
      '<div class="drow"><div class="dk">3&#39; G&gt;A (pos 1)</div><div class="dbarwrap"><div class="dbar" style="width:'+Math.min(100,(d.ga1||0)*100/0.3).toFixed(1)+'%;background:#4f83c2"></div></div><div class="dv">'+(d.ga1==null?'NA':(d.ga1*100).toFixed(1)+'%')+'</div><div class="dp"></div></div>'+
      '<div class="drow"><div class="dk">mean frag len</div><div class="dv" style="flex:1;text-align:left;color:var(--txt2)">'+(d.fraglen==null?'NA':d.fraglen.toFixed(0)+' bp')+'</div></div>';
    }else{dmg='<div class="dsub">aDNA damage</div><div class="nd" style="padding:4px 0">no mapDamage2 output found.</div>';}}
  var fl=s.f.length?'':'<span style="color:#16a34a">no flags ✓</span>';
  var why=s.f.length?'<div class="rsns">'+s.f.map(function(f){var r=flagReason(s,f);return '<span class="rsn '+(FAILF[f]?'bad':'warn')+'" title="'+esc(f)+': '+esc(flagWhy(s,f))+'"><b>'+esc(r[0])+'</b> '+r[1]+'</span>';}).join('')+'</div>':'';
  var annb='';
  (function(){if(s.m.ann_high==null&&s.m.annotated_pct==null)return;
    var warn=s.m.snpeff_warn, tot=s.m.total_variants, badFrac=(warn!=null&&tot)?warn/tot:null;
    var annBad=(s.m.annotated_pct!=null&&s.m.annotated_pct<50)||s.ann_db_error;
    var cls=(annBad||(badFrac!=null&&badFrac>0.5))?' low':'';
    annb='<div class="dsub">Functional annotation (snpEff) <span style="font-weight:400;color:#94a3b8">- impact + annotation QC</span></div>'+
      '<div class="acards">'+
      '<div class="acard'+cls+'"><div class="acard-h"><span class="sname">Impact</span></div>'+
        '<div class="astats">'+
        '<div class="astat"><span class="ak">HIGH</span><span class="av'+(s.m.ann_high>0?' bad':'')+'">'+(s.m.ann_high==null?'NA':Math.round(s.m.ann_high).toLocaleString("en-US"))+'</span></div>'+
        '<div class="astat"><span class="ak">MODERATE</span><span class="av">'+(s.m.ann_moderate==null?"NA":Math.round(s.m.ann_moderate).toLocaleString("en-US"))+'</span></div>'+
        '<div class="astat"><span class="ak">LoF %</span><span class="av">'+(s.m.lof_pct==null?"NA":s.m.lof_pct.toFixed(1)+"%")+'</span></div>'+
        '<div class="astat"><span class="ak">missense/silent</span><span class="av">'+(s.m.missense_silent==null?"NA":s.m.missense_silent.toFixed(2))+'</span></div>'+
        '</div></div>'+
      '<div class="acard'+cls+'"><div class="acard-h"><span class="sname">Annotation QC</span></div>'+
        '<div class="astats">'+
        '<div class="astat"><span class="ak">annotated</span><span class="av'+(annBad?' bad':' good')+'">'+(s.m.annotated_pct==null?"NA":s.m.annotated_pct.toFixed(1)+"%")+'</span></div>'+
        '<div class="astat"><span class="ak">coding</span><span class="av">'+(s.m.coding_pct==null?"NA":s.m.coding_pct.toFixed(1)+"%")+'</span></div>'+
        '<div class="astat"><span class="ak">snpEff warnings</span><span class="av'+(badFrac!=null&&badFrac>0.5?' bad':'')+'">'+(s.m.snpeff_warn==null?"NA":Math.round(s.m.snpeff_warn).toLocaleString("en-US"))+'</span></div>'+
        '</div>'+((annBad||(badFrac!=null&&badFrac>0.5))?'<div class="alow">Low annotation coverage, a snpEff database error, or many snpEff warnings: the reference GFF-to-database build may be wrong for this contig; treat the impact counts and the missense/silent proxy for this sample with suspicion.</div>':'')+
      '</div></div>';})();
  el('modalbody').innerHTML='<div class="dhead"><div><div class="dtitle">'+esc(s.s)+'</div><div class="dmeta" title="'+esc(s.lineage||'')+'">'+esc(s.lineage?linLabel(s.lineage):'lineage NA')+(s.ref?' &middot; mapped to '+esc(s.ref):'')+(s.anc?' &middot; <b style="color:#8a5a12">aDNA</b>':'')+(s.dr?' &middot; DR: '+esc(s.dr):'')+(s.date?' &middot; '+esc(s.date):'')+'</div></div>'+
    '<span class="v '+s.v+'" style="font-size:12px">'+s.v+'</span></div><div class="dflags">'+fl+'</div>'+why+genomeSpark(s)+lin+dmg+annb+
    '<div class="dsub">All metrics <span style="font-weight:400;color:#94a3b8">(bar = position in cohort range · p = percentile)</span></div>'+rows+
    '<div class="dbtns"><button class="btn" id="dexcl"></button></div>';
  var dx=el('dexcl');function setlbl(){dx.textContent=st.excl[s.s]?'✓ excluded - click to keep':'exclude this sample';dx.classList.toggle('prim',!st.excl[s.s]);}
  setlbl();dx.onclick=function(){if(st.excl[s.s])delete st.excl[s.s];else st.excl[s.s]=1;setlbl();renderTable();renderCuration();};
  st._opener=document.activeElement; el('modal').classList.add('open'); setBgInert(true); var mx=el('modalx'); if(mx)mx.focus();}
function closeDetail(){st.detail=null;el('modal').classList.remove('open');setBgInert(false);
  if(st._opener&&st._opener.focus){try{st._opener.focus();}catch(e){}} st._opener=null;}

// ---- genome landscape: samples x reference-position heatmap, multi-track (missing / SNPs / het / indels) ----
