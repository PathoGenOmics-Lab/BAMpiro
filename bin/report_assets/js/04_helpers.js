function esc(s){return String(s).replace(/[&<>"]/g,function(c){return{'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;'}[c];});}
// ---- Icon set (Lucide, MIT): one homogeneous stroke family; currentColor -> auto light/dark ----
var IC={
  printer:'<path d="M6 18H4a2 2 0 0 1-2-2v-5a2 2 0 0 1 2-2h16a2 2 0 0 1 2 2v5a2 2 0 0 1-2 2h-2"/><path d="M6 9V3a1 1 0 0 1 1-1h10a1 1 0 0 1 1 1v6"/><rect x="6" y="14" width="12" height="8" rx="1"/>',
  spark:'<path d="M12 2.5l2.1 5.9 5.9 2.1-5.9 2.1L12 18.5l-2.1-6L4 10.5l5.9-2.1z"/><path d="M19 15l.7 2 2 .7-2 .7-.7 2-.7-2-2-.7 2-.7z"/>',
  grid:'<rect x="3" y="3" width="7" height="7" rx="1"/><rect x="14" y="3" width="7" height="7" rx="1"/><rect x="14" y="14" width="7" height="7" rx="1"/><rect x="3" y="14" width="7" height="7" rx="1"/>',
  sun:'<circle cx="12" cy="12" r="4"/><path d="M12 2v2"/><path d="M12 20v2"/><path d="m4.9 4.9 1.4 1.4"/><path d="m17.7 17.7 1.4 1.4"/><path d="M2 12h2"/><path d="M20 12h2"/><path d="m6.3 17.7-1.4 1.4"/><path d="m19.1 4.9-1.4 1.4"/>',
  moon:'<path d="M12 3a6 6 0 0 0 9 9 9 9 0 1 1-9-9Z"/>',
  panel:'<rect width="18" height="18" x="3" y="3" rx="2"/><path d="M9 3v18"/>',
  sliders:'<line x1="21" x2="14" y1="4" y2="4"/><line x1="10" x2="3" y1="4" y2="4"/><line x1="21" x2="12" y1="12" y2="12"/><line x1="8" x2="3" y1="12" y2="12"/><line x1="21" x2="16" y1="20" y2="20"/><line x1="12" x2="3" y1="20" y2="20"/><line x1="14" x2="14" y1="2" y2="6"/><line x1="8" x2="8" y1="10" y2="14"/><line x1="16" x2="16" y1="18" y2="22"/>',
  maximize:'<polyline points="15 3 21 3 21 9"/><polyline points="9 21 3 21 3 15"/><line x1="21" x2="14" y1="3" y2="10"/><line x1="3" x2="10" y1="21" y2="14"/>',
  minimize:'<polyline points="4 14 10 14 10 20"/><polyline points="20 10 14 10 14 4"/><line x1="14" x2="21" y1="10" y2="3"/><line x1="3" x2="10" y1="21" y2="14"/>',
  basket:'<path d="m15 11-1 9"/><path d="m19 11-4-7"/><path d="M2 11h20"/><path d="m3.5 11 1.6 7.4a2 2 0 0 0 2 1.6h9.8a2 2 0 0 0 2-1.6l1.6-7.4"/><path d="m5 11 4-7"/><path d="m9 11 1 9"/>',
  download:'<path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4"/><polyline points="7 10 12 15 17 10"/><line x1="12" x2="12" y1="15" y2="3"/>',
  search:'<circle cx="11" cy="11" r="8"/><path d="m21 21-4.3-4.3"/>',
  x:'<path d="M18 6 6 18"/><path d="m6 6 12 12"/>',
  chevronDown:'<path d="m6 9 6 6 6-6"/>',
  chevronUp:'<path d="m18 15-6-6-6 6"/>',
  help:'<circle cx="12" cy="12" r="10"/><path d="M9.09 9a3 3 0 0 1 5.83 1c0 2-3 3-3 3"/><path d="M12 17h.01"/>',
  alert:'<path d="m21.73 18-8-14a2 2 0 0 0-3.48 0l-8 14A2 2 0 0 0 4 21h16a2 2 0 0 0 1.73-3Z"/><path d="M12 9v4"/><path d="M12 17h.01"/>',
  ban:'<circle cx="12" cy="12" r="10"/><path d="m4.9 4.9 14.2 14.2"/>',
  github:'<path d="M15 22v-4a4.8 4.8 0 0 0-1-3.5c3 0 6-2 6-5.5.08-1.25-.27-2.48-1-3.5.28-1.15.28-2.35 0-3.5 0 0-1 0-3 1.5-2.64-.5-5.36-.5-8 0C6 2 5 2 5 2c-.3 1.15-.3 2.35 0 3.5A5.4 5.4 0 0 0 4 9c0 3.5 3 5.5 6 5.5-.39.49-.68 1.05-.85 1.65-.17.6-.22 1.23-.15 1.85v4"/><path d="M9 18c-4.51 2-5-2-7-2"/>',
  ext:'<path d="M15 3h6v6"/><path d="M10 14 21 3"/><path d="M18 13v6a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2V8a2 2 0 0 1 2-2h6"/>'
};
function icon(n,cls){return '<svg class="ic'+(cls?' '+cls:'')+'" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round" aria-hidden="true" focusable="false">'+(IC[n]||'')+'</svg>';}
function fillIcons(root){Array.prototype.forEach.call((root||document).querySelectorAll('[data-ic]'),function(e){e.innerHTML=icon(e.getAttribute('data-ic'),e.getAttribute('data-ic-cls')||'');e.removeAttribute('data-ic');});}
function fmt(v,k){if(v==null)return'NA';if(k=='int')return Math.round(v).toLocaleString('en-US');if(k=='pct')return v.toFixed(1);return v.toFixed(2);}
function shortv(v,k){if(v==null)return'';if(k=='int')return Math.round(v).toLocaleString('en-US');return (+v).toPrecision(3);}
function el(id){return document.getElementById(id);}
function visible(){return R.samples.filter(function(s){
  if(st.onlyFlagged&&s.v=='PASS')return false;
  if(st.flagFilter&&s.f.indexOf(st.flagFilter)<0)return false;
  if(st.ancOnly=='anc'&&!s.anc)return false;
  if(st.ancOnly=='mod'&&s.anc)return false;
  if(st.linFilter&&s.lineage!=st.linFilter)return false;
  if(st.q&&s.s.toLowerCase().indexOf(st.q)<0)return false; return true;});}
// Per-column filters for the General Statistics table (table-scoped; do NOT touch the global visible()
// so the plots stay driven by the global filters). A filter string is a numeric operator/range on numeric
// columns (>50, >=50, <10, 5-9, 5..9), otherwise a case-insensitive substring on the displayed cell text.
function colMatchOne(raw,val,txt){
  var q=(raw||'').trim(); if(!q) return true;
  if(typeof val=='number'&&!isNaN(val)){
    var m=q.match(/^(>=|<=|>|<|=)\s*(-?\d+(?:\.\d+)?)$/);
    if(m){var n=parseFloat(m[2]),o=m[1];
      if(o=='>')return val>n; if(o=='>=')return val>=n; if(o=='<')return val<n; if(o=='<=')return val<=n; return val==n;}
    m=q.match(/^(-?\d+(?:\.\d+)?)\s*(?:\.\.|-|to)\s*(-?\d+(?:\.\d+)?)$/);
    if(m){var a=parseFloat(m[1]),b=parseFloat(m[2]); if(a>b){var t=a;a=b;b=t;} return val>=a&&val<=b;}
  }
  return String(txt==null?'':txt).toLowerCase().indexOf(q.toLowerCase())>=0;
}
function colFilterVal(s,k){   // -> [numericValueOrNull, displayText] for column key k
  if(k=='s')return [null,s.s];
  if(k=='v')return [null,s.v];
  if(k=='lineage')return [null,s.lineage||'NA'];
  var v=s.m[k]; return [v,(v==null?'NA':fmt(v,(MET[k]||{}).kind))];
}
function colMatch(s){if(!st.showColF)return true;   // filters apply only while the filter row is shown
  for(var k in st.colf){var raw=st.colf[k]; if(!raw||!raw.trim())continue;
  var pv=colFilterVal(s,k); if(!colMatchOne(raw,pv[0],pv[1]))return false;} return true;}
function colAnyActive(){if(!st.showColF)return false; for(var k in st.colf){if(st.colf[k]&&st.colf[k].trim())return true;} return false;}
// any row filter active (used to surface a "clear filters" affordance so users never lose track of why rows vanished)
function anyFilterActive(){return !!(st.q||st.onlyFlagged||st.flagFilter||st.linFilter||st.ancOnly||colAnyActive());}
function clearAllFilters(){
  st.q=''; st.onlyFlagged=false; st.flagFilter=null; st.linFilter=null; st.ancOnly=null; st.colf={};
  var q=el('q'); if(q)q.value=''; var of=el('of'); if(of)of.checked=false;
  Array.prototype.forEach.call(document.querySelectorAll('#ancseg button'),function(x){x.classList.toggle('on',(x.getAttribute('data-a')||'')=='');});
  renderAll();
}
function dotColor(s){return st.colorBy=='lineage'?linColor(s.lineage):VCOL[s.v];}
// shared colour key for every dot plot (adapts to the QC/lineage colour toggle)
function colorLegend(){
  var sw='display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:5px;vertical-align:-1px;background:';
  var items=(st.colorBy=='lineage'
    ? (R.lineages||[]).map(function(l){return '<span><i style="'+sw+linColor(l)+'"></i>'+esc(l)+'</span>';}).join('')||'<span class="c">no lineage assigned</span>'
    : '<span><i style="'+sw+VCOL.PASS+'"></i>PASS</span><span><i style="'+sw+VCOL.WARN+'"></i>WARN</span><span><i style="'+sw+VCOL.FAIL+'"></i>FAIL</span>');
  return '<div class="legend" style="justify-content:center">'+items+'</div>';
}
// genome landscape: missing-fraction (0 callable -> 1 missing) mapped to a pale->red heat colour
function heatCol(mv){if(mv==null)return TH.cellnull;var a=isDark()?[34,46,60]:[238,244,240],b=[214,64,58];
  return 'rgb('+Math.round(a[0]+(b[0]-a[0])*mv)+','+Math.round(a[1]+(b[1]-a[1])*mv)+','+Math.round(a[2]+(b[2]-a[2])*mv)+')';}
function fmtpos(p){return p>=1e6?(p/1e6).toFixed(2)+' Mb':p>=1e3?Math.round(p/1e3)+' kb':(''+p)+' bp';}
function tip(h,x,y){var t=el('tt'); if(!h){t.style.opacity=0;return;} t.innerHTML=h;
  var w=window.innerWidth; t.style.left=Math.min(x+13,w-t.offsetWidth-10)+'px'; t.style.top=(y+13)+'px'; t.style.opacity=1;}
function setHi(s){st.hi=(st.hi==s?null:s); renderAll();}

function donut(c){var t=(c.PASS+c.WARN+c.FAIL)||1,R0=38,C=2*Math.PI*R0,off=0,segs='';
  segs+='<circle cx="46" cy="46" r="'+R0+'" fill="none" stroke="'+TH.track+'" stroke-width="13"/>';
  [['PASS',VFILL.PASS],['WARN',VFILL.WARN],['FAIL',VFILL.FAIL]].forEach(function(p){var frac=c[p[0]]/t,len=frac*C;
    if(len<=0)return;
    segs+='<circle cx="46" cy="46" r="'+R0+'" fill="none" stroke="'+p[1]+'" stroke-width="13" stroke-dasharray="'+len.toFixed(2)+' '+(C-len).toFixed(2)+'" stroke-dashoffset="'+(-off).toFixed(2)+'" transform="rotate(-90 46 46)"/>'; off+=len;});
  var pct=Math.round(100*c.PASS/t);
  return '<svg width="92" height="92" viewBox="0 0 92 92">'+segs+'<text x="46" y="42" text-anchor="middle" font-size="22" font-weight="700" fill="'+TH.ink+'" letter-spacing="-.5">'+pct+'%</text><text x="46" y="58" text-anchor="middle" font-size="9.5" fill="'+TH.axis+'" letter-spacing=".1em">PASS</text></svg>';}

// ---- Executive summary: cohort KPIs, quality profile and headline findings (for the PI receiving the file)
function _median(vals){var a=vals.filter(function(v){return v!=null;}).sort(function(x,y){return x-y;}); if(!a.length)return null; var m=Math.floor(a.length/2); return a.length%2?a[m]:(a[m-1]+a[m])/2;}
var INS_ICON='<svg viewBox="0 0 24 24" width="14" height="14" fill="currentColor" aria-hidden="true"><path d="M12 2.5l2.1 5.9 5.9 2.1-5.9 2.1L12 18.5l-2.1-6L4 10.5l5.9-2.1z"/></svg>';
// Generic analytical read-out block. title=short caps label; narrHTML=the computed verdict sentence(s);
// chips=optional array of {t,cls,title} rendered as pills (cls: ''|good|warn|bad). Hidden en masse by the header toggle.
function insBox(title, narrHTML, chips){
  var ch=(chips&&chips.length)?('<div class="ins-chips">'+chips.map(function(c){return '<span class="ins-chip'+(c.cls?(' '+c.cls):'')+'"'+(c.title?(' title="'+esc(c.title)+'"'):'')+'>'+c.t+'</span>';}).join('')+'</div>'):'';
  return '<div class="insight"><div class="ins-h"><span class="ins-title">'+INS_ICON+esc(title)+'</span></div><div class="ins-narr">'+narrHTML+'</div>'+ch+'</div>';
}
function _range(vals){var a=vals.filter(function(v){return v!=null;}); return a.length?[Math.min.apply(null,a),Math.max.apply(null,a)]:null;}
/* ---- per-panel analytical read-outs (one insight per panel; each is standalone, computed from R) ---- */
