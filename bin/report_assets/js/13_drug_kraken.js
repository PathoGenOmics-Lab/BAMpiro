function drGColor(gn){ return (gn===1||gn===2)?'#dc2626':(gn===3?'#d97706':((gn===4||gn===5)?'#94a3b8':'#b8c2cf')); }
function drGInk(gn){ return (gn===1||gn===2||gn===3)?'#fff':'#1f2a37'; }   // dark ink on the pale grey (4-5 / unknown) badges so the label stays legible
function drStatus(gns){ var r=false,u=false,n=false; for(var i=0;i<gns.length;i++){var g=gns[i]; if(g===1||g===2)r=true; else if(g===3)u=true; else if(g===4||g===5)n=true;} return r?{t:'R',c:'#dc2626'}:(u?{t:'?',c:'#d97706'}:(n?{t:'&#183;',c:'#94a3b8'}:{t:'',c:'var(--track)'})); }
// ---- Kraken2 taxonomic composition (contamination / host check on the reads before mapping) ----
function renderKraken(){
  var host=el('kraken_body'), sec=el('kraken'); if(!host)return;
  var K=R.kraken;
  if(!(K&&K.samples&&K.samples.length)){ if(sec)sec.style.display='none'; var nv=el('nav-kraken'); if(nv)nv.style.display='none'; return; }
  var rows=K.samples.slice().sort(function(a,b){ return (a.primary?a.primary.pct:0)-(b.primary?b.primary.pct:0); });  // most-contaminated (lowest primary %) first
  // ---- cohort stacked-composition plot: one bar per sample, full taxa breakdown ----
  var agg={}; rows.forEach(function(s){ (s.top||[]).forEach(function(t){ agg[t.name]=(agg[t.name]||0)+t.pct; }); });
  var taxa=Object.keys(agg).sort(function(a,b){ return agg[b]-agg[a]; }).slice(0,8);
  var TAXPAL=['#3b7dd8','#e0544f','#2ea36b','#e0a11f','#8a63c9','#26a0a0','#d06fae','#c98a3b'], OTHERC='#9aa7b6';
  var tcol={}; taxa.forEach(function(t,i){ tcol[t]=TAXPAL[i%TAXPAL.length]; });
  var kW=Math.max(320,(host.clientWidth||760)), kLab=Math.min(150,Math.round(kW*0.28)), kRP=10, kBx=kLab, kBw=kW-kLab-kRP,
      kRowH=Math.max(15,Math.min(22,Math.floor(340/rows.length))), kTop=6, kH=kTop+rows.length*kRowH+20;
  var ksvg='<svg width="'+kW+'" height="'+kH+'" style="display:block;max-width:100%">';
  [0,0.5,1].forEach(function(f){ var x=kBx+f*kBw; ksvg+='<line x1="'+x.toFixed(1)+'" y1="'+kTop+'" x2="'+x.toFixed(1)+'" y2="'+(kTop+rows.length*kRowH).toFixed(1)+'" stroke="'+TH.grid+'"/><text x="'+x.toFixed(1)+'" y="'+(kH-6)+'" text-anchor="'+(f===0?'start':f===1?'end':'middle')+'" font-size="9.5" fill="'+TH.mut+'">'+(f*100)+'%</text>'; });
  rows.forEach(function(s,i){
    var y=kTop+i*kRowH, cx=kBx, used=0, arr=s.top||[];
    ksvg+='<g data-s="'+esc(s.s)+'" style="cursor:pointer"><text x="'+(kLab-6)+'" y="'+(y+kRowH/2+3).toFixed(1)+'" text-anchor="end" font-size="10" fill="'+TH.ink+'">'+esc(s.s.length>20?s.s.slice(0,19)+'…':s.s)+'</text>';
    taxa.forEach(function(t){ var m=null; for(var j=0;j<arr.length;j++){ if(arr[j].name===t){ m=arr[j]; break; } }
      if(m&&m.pct>0){ var w=m.pct/100*kBw; ksvg+='<rect x="'+cx.toFixed(1)+'" y="'+(y+2).toFixed(1)+'" width="'+Math.max(0.4,w).toFixed(1)+'" height="'+(kRowH-4)+'" fill="'+tcol[t]+'"><title>'+esc(s.s)+' · '+esc(t)+' '+m.pct.toFixed(1)+'%</title></rect>'; cx+=w; used+=m.pct; } });
    var other=Math.max(0,(s.classified||0)-used); if(other>0.05){ var wo=other/100*kBw; ksvg+='<rect x="'+cx.toFixed(1)+'" y="'+(y+2).toFixed(1)+'" width="'+wo.toFixed(1)+'" height="'+(kRowH-4)+'" fill="'+OTHERC+'"><title>'+esc(s.s)+' · other classified '+other.toFixed(1)+'%</title></rect>'; cx+=wo; }
    var unc=s.unclassified||0; if(unc>0.05){ var wu=unc/100*kBw; ksvg+='<rect x="'+cx.toFixed(1)+'" y="'+(y+2).toFixed(1)+'" width="'+wu.toFixed(1)+'" height="'+(kRowH-4)+'" fill="'+TH.track+'"><title>'+esc(s.s)+' · unclassified '+unc.toFixed(1)+'%</title></rect>'; }
    ksvg+='</g>';
  });
  ksvg+='</svg>';
  var kleg='<div class="krk-legend">'+taxa.map(function(t){ return '<span><i style="background:'+tcol[t]+'"></i>'+esc(t)+'</span>'; }).join('')+
    '<span><i style="background:'+OTHERC+'"></i>other classified</span><span><i style="background:'+TH.track+'"></i>unclassified</span> <span class="krk-mut">- one bar per sample, worst first; hover a segment for the %. Click a bar/row to highlight the sample everywhere.</span></div>';
  var body=rows.map(function(s){
    var pri=s.primary?s.primary.pct:0, unc=s.unclassified||0, other=Math.max(0,100-pri-unc), lowPri=pri<90, hiUnc=unc>15;
    var bar='<div class="krk-bar">'+
      '<div class="krk-seg" style="width:'+pri.toFixed(1)+'%;background:'+(lowPri?'#e0a11f':'#2ea36b')+'" title="primary: '+esc(s.primary?s.primary.name:'-')+' '+pri.toFixed(1)+'%"></div>'+
      '<div class="krk-seg" style="width:'+other.toFixed(1)+'%;background:#c0704f" title="other classified '+other.toFixed(1)+'%"></div>'+
      '<div class="krk-seg" style="width:'+unc.toFixed(1)+'%;background:var(--track)" title="unclassified '+unc.toFixed(1)+'%"></div></div>';
    return '<tr class="hit" data-s="'+esc(s.s)+'"'+(st.hi==s.s?' style="background:'+TH.hl+'"':'')+'><td class="s">'+esc(s.s)+'</td>'+
      '<td style="text-align:left"><b>'+esc(s.primary?s.primary.name:'-')+'</b></td>'+
      '<td'+(lowPri?' style="color:var(--fail);font-weight:600"':'')+'>'+pri.toFixed(1)+'</td>'+
      '<td style="text-align:left">'+((s.secondary&&s.secondary.pct>=1)?esc(s.secondary.name)+' <span class="krk-mut">'+s.secondary.pct.toFixed(1)+'%</span>':'<span class="krk-mut">-</span>')+'</td>'+
      '<td'+(hiUnc?' style="color:var(--warn)"':'')+'>'+unc.toFixed(1)+'</td>'+
      '<td class="krk-barcell">'+bar+'</td></tr>';
  }).join('');
  host.innerHTML='<div class="krk-chart">'+kleg+'<div class="krk-plotscroll">'+ksvg+'</div></div>'+
    '<div class="gtable" style="margin-top:14px"><table class="krktable"><thead><tr><th class="s">Sample</th><th style="text-align:left">Primary taxon</th><th>Primary %</th><th style="text-align:left">Top other</th><th>Unclass. %</th><th style="text-align:left">Composition</th></tr></thead><tbody>'+body+'</tbody></table></div>';
  Array.prototype.forEach.call(host.querySelectorAll('[data-s]'),function(e){e.onclick=function(){setHi(e.getAttribute('data-s'));};});
}
var drState={q:''};
function renderDrug(){
  var host=el('drug_body'), sec=el('drug'); if(!host)return;
  var D=R.dr;
  if(!(D&&D.calls&&D.calls.length)){ if(sec)sec.style.display='none'; var nv=el('nav-drug'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  var samples=D.samples, drugs=D.drugs, calls=D.calls;
  var cell={}, cmut={};
  calls.forEach(function(c){ (cell[c.s]=cell[c.s]||{}); (cell[c.s][c.drug]=cell[c.s][c.drug]||[]).push(c.gn); (cmut[c.s]=cmut[c.s]||{}); (cmut[c.s][c.drug]=cmut[c.s][c.drug]||[]).push(c); });
  var mx='<div class="dr-mxwrap"><table class="drmx"><thead><tr><th class="dr-corner">sample \\ drug</th>'+
    drugs.map(function(dr){return '<th class="dr-hcell" title="'+esc(dr)+'"><span class="dr-h">'+esc(dr)+'</span></th>';}).join('')+'</tr></thead><tbody>'+
    samples.map(function(s){ return '<tr><th class="dr-row" title="'+esc(s)+'">'+esc(s)+'</th>'+drugs.map(function(dr){
      var gns=(cell[s]||{})[dr];
      if(!gns) return '<td class="drmx-cell" title="'+esc(s)+' &#183; '+esc(dr)+': no mutation detected"></td>';
      var st=drStatus(gns);
      var muts=((cmut[s]||{})[dr]||[]).map(function(c){return c.gene+' '+c.mutation+(c.gn?(' (WHO '+c.gn+')'):'');}).join('; ');
      return '<td class="drmx-cell" style="background:'+st.c+'" title="'+esc(s)+' &#183; '+esc(dr)+' &#8212; '+esc(muts)+'"><b>'+st.t+'</b></td>';
    }).join('')+'</tr>'; }).join('')+'</tbody></table></div>';
  host.innerHTML=
    '<div class="dr-legend">'+
      '<span title="WHO groups 1-2: associated with resistance"><i style="background:#dc2626"></i>1&#8211;2 associated with R</span>'+
      '<span title="WHO group 3: uncertain significance"><i style="background:#d97706"></i>3 uncertain</span>'+
      '<span title="WHO groups 4-5: not associated with resistance"><i style="background:#94a3b8"></i>4&#8211;5 not associated</span>'+
      '<span class="c">Cell = worst grade per drug (R / ? / &#183;). A genomic screen (pathotypr, WHO catalogue, H37Rv numbering), not a clinical DST result.</span></div>'+
    mx+
    '<div class="dr-controls"><input id="drq" class="dyn-search" type="search" placeholder="filter by sample / drug / gene / mutation..." value="'+esc(drState.q)+'"><button class="dyn-btn" id="drdl" title="Download every resistance call as a TSV">'+icon('download')+'download calls (TSV)</button><span class="dyn-count" id="drcount"></span></div>'+
    '<div class="epitbl-wrap"><table class="epitbl" id="drtable"></table></div>';
  function draw(){
    var q=drState.q.toLowerCase();
    var rows=calls.filter(function(c){ return !q||(c.s.toLowerCase().indexOf(q)>=0)||(c.drug.toLowerCase().indexOf(q)>=0)||(c.gene.toLowerCase().indexOf(q)>=0)||(c.mutation.toLowerCase().indexOf(q)>=0); });
    rows=rows.slice().sort(function(a,b){ if(a.s!==b.s) return a.s<b.s?-1:1; return (a.gn||9)-(b.gn||9); });
    el('drcount').textContent=rows.length+' call(s)'+(rows.length>600?' · showing first 600 (download for all)':'');
    var h='<thead><tr><th>Sample</th><th>Drug</th><th>Gene</th><th>Mutation (H37Rv)</th><th>WHO grade</th><th>AF</th><th>DP</th></tr></thead><tbody>';
    if(!rows.length) h+='<tr><td colspan="7" class="c" style="padding:18px;text-align:center">no call matches the filter.</td></tr>';
    h+=rows.slice(0,600).map(function(c){
      return '<tr><td>'+esc(c.s)+'</td><td><b>'+esc(c.drug)+'</b></td><td>'+esc(c.gene)+geneRvTag(c.gene)+'</td><td class="epitbl-r">'+esc(c.mutation)+'</td>'+
        '<td><span class="dr-badge" style="background:'+drGColor(c.gn)+';color:'+drGInk(c.gn)+'" title="'+esc(c.marker||'')+'">'+esc(c.grade||'?')+'</span></td>'+
        '<td>'+(c.af==null?'':c.af.toFixed(2))+'</td><td>'+(c.dp==null?'':c.dp)+'</td></tr>';
    }).join('')+'</tbody>';
    el('drtable').innerHTML=h;
  }
  el('drq').oninput=function(){ drState.q=this.value; draw(); };
  el('drdl').onclick=function(){
    var hdr=['sample','drug','gene','mutation_h37rv','who_grade','marker','af','dp'];
    var lines=[hdr.join('\t')];
    calls.forEach(function(c){ lines.push([c.s,c.drug,c.gene,c.mutation,c.grade,c.marker,(c.af==null?'':c.af),(c.dp==null?'':c.dp)].join('\t')); });
    dl(lines.join('\n')+'\n','drug_resistance.tsv','text/tab-separated-values');
  };
  draw();
}

// ---- Dose x treatment: per-group dose distribution + Kruskal-Wallis test (whole cohort) ----
function renderDoseTx(){
  var host=el('dosetx_body'), sec=el('dosetx'), cap=el('dosetx_caption'), nv=el('nav-dosetx');
  if(!host)return;
  function hide(){ if(sec)sec.style.display='none'; if(nv)nv.style.display='none'; }
  var meta=R.sample_meta, txf=meta&&meta.tx_field;
  if(!(MET.dose&&txf)){ hide(); return; }                        // need a dose metric AND a treatment column
  var groups={}, order=[];
  R.samples.forEach(function(s){                                 // whole cohort (independent of the live filter)
    var d=s.m.dose; if(d==null)return;
    var tx=(meta.rows[s.s]||{})[txf]; if(!tx||tx==='NA'||tx==='.'||tx==='-')return;
    if(!Object.prototype.hasOwnProperty.call(groups,tx)){groups[tx]=[];order.push(tx);}
    groups[tx].push({s:s.s,d:d,v:s.v});
  });
  var cats=order.filter(function(c){return groups[c].length>=1;});
  var withData=cats.filter(function(c){return groups[c].length>=2;});
  if(cats.length<2||withData.length<2){ hide(); return; }        // need >=2 groups with enough dosed samples
  if(sec)sec.style.display=''; if(nv)nv.style.display='';
  var allD=[]; cats.forEach(function(c){groups[c].forEach(function(o){allD.push(o.d);});});
  var lo=Math.min.apply(null,allD), hi=Math.max.apply(null,allD); if(hi<=lo)hi=lo+1;
  var W=Math.max(360,(host.clientWidth||760)), labW=Math.min(160,Math.round(W*0.26)), rp=14,
      plotW=W-labW-rp, rowH=54, top=10, H=top+cats.length*rowH+30;
  function X(d){return labW+(d-lo)/(hi-lo)*plotW;}
  var pal=['#0e8ba8','#a855c9','#e0a11f','#2ea36b','#e0544f','#6b7280'];
  function srt(a){return a.slice().sort(function(x,y){return x-y;});}
  function med(a){var so=srt(a),n=so.length;return n%2?so[(n-1)/2]:(so[n/2-1]+so[n/2])/2;}
  function q(a,p){var so=srt(a),i=(so.length-1)*p,l=Math.floor(i),h=Math.ceil(i);return l==h?so[l]:so[l]+(so[h]-so[l])*(i-l);}
  var svg='<svg width="'+W+'" height="'+H+'" style="display:block;max-width:100%">';
  [0,0.25,0.5,0.75,1].forEach(function(f){var x=labW+f*plotW; svg+='<line x1="'+x.toFixed(1)+'" y1="'+top+'" x2="'+x.toFixed(1)+'" y2="'+(top+cats.length*rowH).toFixed(1)+'" stroke="'+TH.grid+'"/><text x="'+x.toFixed(1)+'" y="'+(H-14)+'" text-anchor="middle" font-size="9.5" fill="'+TH.mut+'">'+shortv(lo+f*(hi-lo),'float')+'</text>';});
  svg+='<text x="'+(labW+plotW/2)+'" y="'+(H-2)+'" text-anchor="middle" font-size="10.5" fill="'+TH.mut+'">'+esc(MET.dose.label)+'</text>';
  cats.forEach(function(c,i){
    var arr=groups[c].map(function(o){return o.d;}), cy=top+i*rowH+rowH/2, col=pal[i%pal.length];
    var gactive=(txf&&st.metaFilter&&st.metaFilter[txf]===c);   // the group label toggles a cohort filter by that treatment
    svg+='<rect class="dtx-glabel" data-g="'+esc(c)+'" x="0" y="'+(top+i*rowH)+'" width="'+labW+'" height="'+rowH+'" rx="4" fill="#0e8ba8" fill-opacity="'+(gactive?0.12:0)+'" style="cursor:pointer"><title>click to filter the whole report by '+esc(c)+'</title></rect>';
    if(arr.length>=2){var q1=q(arr,0.25),q3=q(arr,0.75);
      svg+='<rect x="'+X(q1).toFixed(1)+'" y="'+(cy-11)+'" width="'+Math.max(1,X(q3)-X(q1)).toFixed(1)+'" height="22" rx="3" fill="'+col+'" fill-opacity="0.14" stroke="'+col+'" stroke-opacity="0.5"/>';}
    var m=med(arr); svg+='<line x1="'+X(m).toFixed(1)+'" y1="'+(cy-12)+'" x2="'+X(m).toFixed(1)+'" y2="'+(cy+12)+'" stroke="'+col+'" stroke-width="2"/>';
    var seen={};
    groups[c].forEach(function(o){var x=X(o.d),key=Math.round(x/6),off=((seen[key]=(seen[key]||0)+1)-1),dy=((off%2)?1:-1)*Math.ceil(off/2)*5;
      dy=Math.max(-20,Math.min(20,dy));   // keep a crowded x-bin (many samples at one dose) inside its own row band
      var hi=(st.hi==o.s);                 // clickable: highlight the sample everywhere (setHi)
      svg+='<circle class="dtx-dot" data-s="'+esc(o.s)+'" cx="'+x.toFixed(1)+'" cy="'+(cy+dy).toFixed(1)+'" r="'+(hi?5.4:3.4)+'" fill="'+(o.v=='FAIL'?'#e0544f':col)+'" fill-opacity="0.92" stroke="'+(hi?TH.ink:'#fff')+'" stroke-width="'+(hi?1.6:0.6)+'" style="cursor:pointer"><title>'+esc(o.s)+' · '+esc(c)+' · dose '+o.d+' \u00b7 click to highlight</title></circle>';});
    svg+='<text x="'+(labW-8)+'" y="'+(cy-1)+'" text-anchor="end" font-size="11" font-weight="'+(gactive?'700':'400')+'" fill="'+(gactive?'#0e8ba8':TH.ink)+'" style="pointer-events:none">'+esc(c.length>20?c.slice(0,19)+'…':c)+'</text>'+
      '<text x="'+(labW-8)+'" y="'+(cy+12)+'" text-anchor="end" font-size="9.5" fill="'+TH.mut+'" style="pointer-events:none">n='+arr.length+' · med '+shortv(m,'float')+'</text>';
  });
  svg+='</svg>';
  host.innerHTML=svg;
  Array.prototype.forEach.call(host.querySelectorAll('.dtx-dot'),function(d){d.onclick=function(){setHi(d.getAttribute('data-s'));};});
  Array.prototype.forEach.call(host.querySelectorAll('.dtx-glabel'),function(r){r.onclick=function(){var g=r.getAttribute('data-g'); applyMetaFilter(txf,(st.metaFilter&&st.metaFilter[txf]===g)?'':g);};});
  var kw=kruskalWallis(withData.map(function(c){return groups[c].map(function(o){return o.d;});}));
  if(cap){
    if(kw){var sig=kw.p<0.05;
      cap.innerHTML='<b>Kruskal–Wallis</b> H = '+kw.H.toFixed(2)+' · p = '+pfmt(kw.p)+' · '+kw.k+' groups, N = '+kw.N+
        ' \u00b7 '+(sig?'<span class="sc-sig">dose differs across treatment groups</span>':'no significant difference in dose across groups')+
        '. <span class="krk-mut">Non-parametric rank test over the full cohort; groups with &lt; 2 dosed samples are drawn but not tested. Click a group name to filter the whole report by it; click a point to highlight that sample.</span>';
    } else cap.innerHTML='<span class="krk-mut">Not enough dosed samples per group to test.</span>';
  }
}

// ---- Variant x dose: per-variant allele frequency correlated with dose (Spearman + BH-FDR) ----
var vardoseSelKey=null, vardoseSort={k:'ar',asc:false}, _vdCache, _vdDone=false;
function _vkey(t){return t.r.pos+'|'+t.r.alt;}
function _vdCompute(){   // cohort-fixed scan -> memoized: identical on every renderAll, so compute once; only the DOM render below is cheap
  if(_vdDone)return _vdCache; _vdDone=true; _vdCache=null;
  var M=R.snp_matrix;
  if(!(MET.dose&&M&&M.rows&&M.rows.length&&M.samples&&M.samples.length))return null;
  var dbys={}; R.samples.forEach(function(s){ if(s.m.dose!=null)dbys[s.s]=s.m.dose; });
  var didx=[],dose=[]; M.samples.forEach(function(sid,i){ if(dbys[sid]!=null){didx.push(i);dose.push(dbys[sid]);} });
  if(didx.length<6)return null;                                 // need enough dosed samples to correlate
  var tests=[];
  M.rows.forEach(function(r){                                   // AF absent -> 0 (reference); need >= 3 carriers
    var af=[],carriers=0; for(var k=0;k<didx.length;k++){var c=r.cells[didx[k]],a=c?c[0]:0; af.push(a); if(a>0)carriers++;}
    if(carriers<3)return;
    var rho=spearman(af,dose); if(rho==null)return;             // null when the AF vector has no variance (fixed site)
    tests.push({r:r,af:af,rho:rho,p:spearmanP(rho,didx.length),n:didx.length,carriers:carriers});
  });
  if(!tests.length)return null;
  var qs=bhFDR(tests.map(function(t){return t.p;})); tests.forEach(function(t,i){t.q=qs[i];});
  tests.sort(function(a,b){ return (Math.abs(b.rho)-Math.abs(a.rho))||(a.q-b.q); });
  var nsig=0; tests.forEach(function(t){ if(t.q<=0.05)nsig++; });
  _vdCache={M:M,didx:didx,dose:dose,tests:tests,nsig:nsig};
  return _vdCache;
}
function renderVarDose(){
  var host=el('vardose_body'), sec=el('vardose'), cap=el('vardose_caption'), nv=el('nav-vardose');
  if(!host)return;
  var D=_vdCompute();
  if(!D){ if(sec)sec.style.display='none'; if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display=''; if(nv)nv.style.display='';
  var M=D.M, didx=D.didx, dose=D.dose, tests=D.tests, nsig=D.nsig;
  var vq=(st.vardoseq||'').toLowerCase();   // search box (gene / position / amino acid), filters the view only
  var flt=vq?tests.filter(function(t){var r=t.r;
    return (r.gene&&r.gene.toLowerCase().indexOf(vq)>=0)||String(r.pos).indexOf(vq)>=0||(r.aa&&r.aa.toLowerCase().indexOf(vq)>=0)||(r.pos_h37rv&&(''+r.pos_h37rv).toLowerCase().indexOf(vq)>=0);}):tests;
  // display order: sortable by any column; default is |rho| desc (the memoised base order)
  function keyval(t){var k=vardoseSort.k; return k=='carriers'?t.carriers:(k=='rho'?t.rho:(k=='p'?t.p:(k=='q'?t.q:Math.abs(t.rho))));}
  var ord=flt.slice().sort(function(a,b){var d=keyval(a)-keyval(b); if(!d)d=Math.abs(b.rho)-Math.abs(a.rho); return vardoseSort.asc?d:-d;});
  if(!ord.length){ host.innerHTML='<div class="pad nd">no variant matches &#8220;'+esc(st.vardoseq)+'&#8221;</div>';
    if(cap)cap.innerHTML='<span class="krk-mut">No variant matches the search. '+tests.length+' variant(s) were tested against dose; clear the box to see them.</span>'; return; }
  // selected variant tracked by a stable key so it survives re-sorts / filtering
  var sel=null,selIdx=-1;
  for(var si=0;si<ord.length;si++){ if(_vkey(ord[si])===vardoseSelKey){ sel=ord[si]; selIdx=si; break; } }
  if(!sel){ sel=ord[0]; selIdx=0; vardoseSelKey=_vkey(sel); }
  // scatter of the selected variant: dose (x) vs AF (y); points click through to highlight the sample
  var W=Math.min(440,Math.max(300,Math.round((host.clientWidth||760)*0.42))),Hs=250,padL=44,padB=32,padT=12,padR=12;
  var dlo=Math.min.apply(null,dose),dhi=Math.max.apply(null,dose); if(dhi<=dlo)dhi=dlo+1;
  function X(d){return padL+(d-dlo)/(dhi-dlo)*(W-padL-padR);}
  function Y(a){return Hs-padB-a*(Hs-padB-padT);}
  var svg='<svg width="'+W+'" height="'+Hs+'" style="display:block;max-width:100%">';
  [0,0.25,0.5,0.75,1].forEach(function(t){var y=Y(t);svg+='<line x1="'+padL+'" y1="'+y.toFixed(1)+'" x2="'+(W-padR)+'" y2="'+y.toFixed(1)+'" stroke="'+TH.grid+'"/><text x="'+(padL-6)+'" y="'+(y+3).toFixed(1)+'" text-anchor="end" font-size="9" fill="'+TH.mut+'">'+t.toFixed(2).replace(/^0/,'')+'</text>';});
  [0,0.5,1].forEach(function(t){var d=dlo+t*(dhi-dlo),x=X(d);svg+='<text x="'+x.toFixed(1)+'" y="'+(Hs-padB+13)+'" text-anchor="middle" font-size="9" fill="'+TH.mut+'">'+shortv(d,'float')+'</text>';});
  for(var k=0;k<didx.length;k++){var sid=M.samples[didx[k]],hi=(st.hi==sid),x=X(dose[k]),y=Y(sel.af[k]);
    svg+='<circle class="vd-dot" data-s="'+esc(sid)+'" cx="'+x.toFixed(1)+'" cy="'+y.toFixed(1)+'" r="'+(hi?5.4:3.6)+'" fill="#0e8ba8" fill-opacity="0.82" stroke="'+(hi?TH.ink:'#fff')+'" stroke-width="'+(hi?1.6:0.6)+'" style="cursor:pointer"><title>'+esc(sid)+' · dose '+shortv(dose[k],'float')+' · AF '+sel.af[k].toFixed(2)+' \u00b7 click to highlight</title></circle>';}
  svg+='<text x="'+((padL+W-padR)/2).toFixed(1)+'" y="'+(Hs-3)+'" text-anchor="middle" font-size="10" fill="'+TH.mut+'">'+esc(MET.dose.label)+'</text>'+
    '<text x="11" y="'+((padT+Hs-padB)/2).toFixed(1)+'" text-anchor="middle" font-size="10" fill="'+TH.mut+'" transform="rotate(-90 11 '+((padT+Hs-padB)/2).toFixed(1)+')">allele frequency</text></svg>';
  function vlabel(r){return '<b>'+esc(r.gene||r.contig)+'</b>'+geneRvTag(r.gene)+' '+refPos(r.pos,r.pos_h37rv)+' '+esc(r.ref)+'&#8594;'+esc(r.alt);}
  function rcol(rho){return rho>0?'#c0453b':'#3f6bbf';}
  var MAXT=200, shown=ord.slice(0,MAXT);
  var rowsH=shown.map(function(t){var k=_vkey(t);
    return '<tr class="vd-row'+(k===vardoseSelKey?' on':'')+'" data-k="'+esc(k)+'"><td class="vd-lbl">'+vlabel(t.r)+(t.r.aa?(' <span class="snpmx-aa">'+aaDual(t.r.aa,t.r.aa_h37rv)+'</span>'):'')+'</td>'+
      '<td class="vd-num">'+t.carriers+'/'+t.n+'</td>'+
      '<td class="vd-num" style="color:'+rcol(t.rho)+'">'+(t.rho>0?'+':'&#8722;')+Math.abs(t.rho).toFixed(2)+'</td>'+
      '<td class="vd-num">'+pfmt(t.p)+'</td>'+
      '<td class="vd-num"'+(t.q<=0.05?' style="font-weight:600"':'')+'>'+pfmt(t.q)+(t.q<=0.05?' <span class="sc-sig">*</span>':'')+'</td></tr>';
  }).join('');
  function th(sk,lbl,cls){var on=(vardoseSort.k==sk),ar=on?(vardoseSort.asc?' &#9650;':' &#9660;'):''; return '<th class="vd-sortable'+(cls?(' '+cls):'')+'" data-sk="'+sk+'"'+(on?' style="color:var(--accent)"':'')+'>'+lbl+ar+'</th>';}
  host.innerHTML='<div class="vd-wrap"><div class="vd-scatter">'+
    '<div class="vd-selhdr">'+vlabel(sel.r)+(sel.r.aa?(' <span class="snpmx-aa">'+aaDual(sel.r.aa,sel.r.aa_h37rv)+'</span>'):'')+
    '<div class="vd-selstat">Spearman &#961; = <b style="color:'+rcol(sel.rho)+'">'+(sel.rho>0?'+':'&#8722;')+Math.abs(sel.rho).toFixed(2)+'</b> · p = '+pfmt(sel.p)+' · q = '+pfmt(sel.q)+' · '+sel.carriers+'/'+sel.n+' carriers</div></div>'+svg+'</div>'+
    '<div class="vd-tablewrap"><table class="vd-table"><thead><tr><th class="vd-lbl">Variant</th>'+th('carriers','carriers','vd-num')+th('rho','&#961;','vd-num')+th('p','p','vd-num')+th('q','q (FDR)','vd-num')+'</tr></thead><tbody>'+rowsH+'</tbody></table>'+
    (ord.length>MAXT?('<div class="krk-mut" style="padding:6px 4px">showing the top '+MAXT+' of '+ord.length+' tested variants</div>'):'')+'</div></div>';
  Array.prototype.forEach.call(host.querySelectorAll('.vd-row'),function(tr){tr.onclick=function(){vardoseSelKey=tr.getAttribute('data-k');renderVarDose();};});
  Array.prototype.forEach.call(host.querySelectorAll('.vd-dot'),function(d){d.onclick=function(){setHi(d.getAttribute('data-s'));};});
  Array.prototype.forEach.call(host.querySelectorAll('.vd-sortable'),function(h){h.onclick=function(){var sk=h.getAttribute('data-sk');
    if(vardoseSort.k==sk)vardoseSort.asc=!vardoseSort.asc; else{vardoseSort.k=sk; vardoseSort.asc=(sk=='p'||sk=='q');}  // p/q default ascending (most significant first)
    renderVarDose();};});
  if(cap)cap.innerHTML=tests.length+' variant(s) tested against dose (&#8805; 3 carriers) · <b'+(nsig?' class="sc-sig"':'')+'>'+nsig+' significant at FDR q &#8804; 0.05</b>'+(vq?(' · showing <b>'+ord.length+'</b> matching &#8220;'+esc(st.vardoseq)+'&#8221;'):'')+'. <span class="krk-mut">Spearman rank correlation of per-sample allele frequency (0 where the site is reference) vs dose over the '+didx.length+' dosed samples; Benjamini–Hochberg q across all tested variants. Click a row to plot it, a header to sort, a point to highlight the sample. Complements the Dose × treatment test.</span>';
}

// ---- Gene conversion: candidate tracts, and the evidence that decides whether to believe them ----
// A tract is a run of diagnostic sites carrying the DONOR paralog's alleles. Two other things produce
// exactly that picture: an ordinary substitution that happens to match the paralog, and a read that
// arrived from the donor in the first place. bin/gconv_model.py weighs a tract against both and reports
// a log10 Bayes factor, so the panel leads with that and puts the observable evidence next to it: how
// fixed the donor allele is inside the tract (donor_af_in), whether it also turns up outside (
// donor_af_outside), and whether one molecule carries donor alleles on one side of a breakpoint and
// acceptor alleles on the other, in cis (breakpoint_reads) - the piece a mismapping cannot fake.
var gconvState={q:'',v:'',sk:'bf',asc:false,oneper:true};
var GCONV_MAXROWS=400;
var GCONV_V=[{k:'gene_conversion',lab:'gene conversion',c:'#2ea36b',r:0,
              tip:'a tract explains the reads far better than an independent substitution or reads arriving from the donor'},
             {k:'ambiguous',lab:'ambiguous',c:'#e0a11f',r:1,
              tip:'reported but not called; the reason says what came closest. Often a tract covering every diagnostic site of its locus, which has the same likelihood as every read having come from the donor'},
             {k:'mismapping',lab:'mismapping',c:'#e0544f',r:2,
              tip:'the locus is explained by a fitted fraction of reads arriving from the donor, with nothing left for a tract to account for'},
             {k:'coverage_shift',lab:'coverage shift',c:'#8b6fd6',r:3,
              tip:'the acceptor lost its reads to the donor over a run of sites. Consistent with a conversion longer than the library insert, and equally with a deletion. Not a conversion call'},
             {k:'reference_artifact',lab:'reference artifact',c:'#7a8794',r:4,
              tip:'present in nearly every sample of the cohort. The reference being wrong here, or the aligner doing this to everybody, explains that more simply than the same conversion arising in every isolate. In a CLONAL cohort it may instead be shared ancestry, which recurrence alone cannot distinguish. Only a cohort can make this call at all'},
             {k:'reciprocal_exchange',lab:'reciprocal exchange',c:'#c77d3a',r:6,
              tip:'the donor carries the ACCEPTOR\'s bases over the same stretch, so both copies changed. That is an exchange between them rather than one being overwritten, and gene conversion is non-reciprocal by definition'},
             {k:'reference_derived',lab:'reference derived',c:'#4a90b8',r:5,
              tip:'an outgroup says the REFERENCE carries the derived base over this stretch and the reads carry the ancestral one. The sample changed nothing; the finding belongs to the reference. Without an outgroup this is the same picture as a conversion'}];
// The settings that produced these verdicts, shown beside them. A panel that displays a verdict
// without saying under which rules cannot be checked against another run, and the global run
// header only carries the reference and the container.
function gconvSettings(){
  var p=R.provenance||{}, keys=[], out=[];
  for(var k in p){ if(k.indexOf('gconv_')===0)keys.push(k); }
  if(!keys.length)return '';
  keys.sort();
  for(var i=0;i<keys.length;i++)out.push(esc(keys[i].replace('gconv_',''))+' '+esc(p[keys[i]]));
  return '<div class="gcv-settings" title="Every setting here moves the verdicts above. The output TSV carries the full list as # header lines.">settings &#183; '+out.join(' &#183; ')+'</div>';
}
function gconvDef(v){for(var i=0;i<GCONV_V.length;i++){if(GCONV_V[i].k===v)return GCONV_V[i];}
  return {k:v||'',lab:(v||'unknown').replace(/_/g,' '),c:'#94a3b8',r:3,tip:'verdict not recognised'};}
function gconvId(t){return t.s+'|'+t.pair+'|'+t.start;}
function gconvLocus(t){return esc(t.contig||'?')+' <span class="gcv-from">&#8592;</span> '+esc(t.donor||'?');}
function gconvNum(v,dp){return v==null?'<span class="gcv-na">n/a</span>':(+v).toFixed(dp);}
// donor allele fraction as a bar + the number, so "how fixed" is readable at a glance per row
function gconvAF(v,col,why){
  if(v==null)return '<span class="gcv-na" title="'+esc(why)+'">n/a</span>';
  var w=Math.max(0,Math.min(1,v))*100;
  return '<span class="gcv-ev"><span class="gcv-bar"><span style="width:'+w.toFixed(0)+'%;background:'+col+'"></span></span>'+v.toFixed(2)+'</span>';
}
function renderGconv(){
  var host=el('gconv_body'), sec=el('gconv'), cap=el('gconv_caption'), nv=el('nav-gconv'), sb=el('gconvq');
  if(!host)return;
  var G=R.gconv;
  if(!(G&&G.tracts&&G.tracts.length)){ if(sec)sec.style.display='none'; if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display=''; if(nv)nv.style.display='';
  var vis={}; visible().forEach(function(s){vis[s.s]=1;});   // the cohort filter drives this panel too
  var rows=G.tracts.filter(function(t){return vis[t.s];});
  if(sb){ if(sb.value!==gconvState.q)sb.value=gconvState.q;
    sb.oninput=function(){ gconvState.q=this.value.trim(); draw(); }; }
  var hasRep=rows.some(function(t){return t.rep!=null;});
  function match(t){
    // A gene family reports one converted stretch once per relationship. Showing them all makes
    // one event look like three findings, so by default only the row that stands for the event
    // and its source is listed, and the rest are one click away.
    if(hasRep&&gconvState.oneper&&t.rep!==1)return false;
    if(gconvState.v&&t.verdict!==gconvState.v)return false;
    var q=gconvState.q.toLowerCase(); if(!q)return true;
    return (t.s.toLowerCase().indexOf(q)>=0)||((t.contig||'').toLowerCase().indexOf(q)>=0)||
      ((t.donor||'').toLowerCase().indexOf(q)>=0)||((t.verdict||'').replace(/_/g,' ').indexOf(q)>=0)||
      ((t.reason||'').toLowerCase().indexOf(q)>=0)||((''+t.pair).indexOf(q)>=0);
  }
  function keyof(t,k){
    if(k=='s')return t.s;
    if(k=='locus')return (t.contig||'')+'|'+(t.donor||'');
    if(k=='verdict')return gconvDef(t.verdict).r;
    var v=t[k]; return v==null?-1:v;   // an unmeasurable value sorts to one end, never silently as 0
  }
  function draw(){
    var sel=rows.filter(match), selSet={}; sel.forEach(function(t){selSet[gconvId(t)]=1;});
    var counts={}, nbp=0; rows.forEach(function(t){counts[t.verdict]=(counts[t.verdict]||0)+1; if(t.bp_reads>0)nbp++;});
    // ---- verdict chips: the judgement first, and a click filters the table to it ----
    var chips=GCONV_V.map(function(d){ return '<span class="gcv-vchip'+(gconvState.v===d.k?' on':'')+'" data-v="'+d.k+
      '" title="'+esc(d.tip)+'"><i style="background:'+d.c+'"></i>'+d.lab+' <b>'+(counts[d.k]||0)+'</b></span>'; }).join('');
    // ---- evidence map: bounded-ness (x) vs fixedness (y), ringed when a read crosses a breakpoint ----
    // Descriptive, not decisive: these are the two fractions a person can go and check in the BAM. The
    // model works off the reads and their qualities, so no line here is a decision boundary and none
    // is drawn. What the model concluded is the colour of the point and the BF column below.
    var W=Math.max(340,Math.min(820,(host.clientWidth||760))), H=250, padL=48, padR=92, padT=14, padB=36;
    var plotW=W-padL-padR, plotH=H-padT-padB, laneX=padL+plotW+44;
    function X(v){return padL+Math.max(0,Math.min(1,v))*plotW;}
    function Y(v){return padT+(1-Math.max(0,Math.min(1,v)))*plotH;}
    function ax(f){return f.toFixed(2).replace(/^0/,'');}
    var svg='<svg width="'+W+'" height="'+H+'" style="display:block;max-width:100%">';
    [0,0.25,0.5,0.75,1].forEach(function(f){
      svg+='<line x1="'+X(f).toFixed(1)+'" y1="'+padT+'" x2="'+X(f).toFixed(1)+'" y2="'+(padT+plotH)+'" stroke="'+TH.grid+'"/>'+
        '<text x="'+X(f).toFixed(1)+'" y="'+(padT+plotH+13)+'" text-anchor="middle" font-size="9" fill="'+TH.mut+'">'+ax(f)+'</text>'+
        '<line x1="'+padL+'" y1="'+Y(f).toFixed(1)+'" x2="'+(padL+plotW)+'" y2="'+Y(f).toFixed(1)+'" stroke="'+TH.grid+'"/>'+
        '<text x="'+(padL-6)+'" y="'+(Y(f)+3).toFixed(1)+'" text-anchor="end" font-size="9" fill="'+TH.mut+'">'+ax(f)+'</text>';});
    svg+='<line x1="'+(laneX-20)+'" y1="'+padT+'" x2="'+(laneX-20)+'" y2="'+(padT+plotH)+'" stroke="'+TH.grid+'" stroke-dasharray="2 3"/>'+
      '<text x="'+laneX+'" y="'+(padT+plotH+13)+'" text-anchor="middle" font-size="9" fill="'+TH.mut+'">n/a</text>';
    rows.forEach(function(t,i){
      if(t.af_in==null)return;                                     // nothing to place on the fixedness axis
      var d=gconvDef(t.verdict), on=selSet[gconvId(t)], hot=(st.hi==t.s), bp=(t.bp_reads>0);
      var j=((i*2654435761)%997)/997-0.5;                          // deterministic jitter: tracts pile up on identical values
      var cx=(t.af_out==null?laneX:X(t.af_out))+j*7, cy=Y(t.af_in)+j*7;
      svg+='<circle class="gcv-dot" data-s="'+esc(t.s)+'" cx="'+cx.toFixed(1)+'" cy="'+cy.toFixed(1)+'" r="'+(hot?6.2:(bp?4.9:3.6))+
        '" fill="'+d.c+'" fill-opacity="'+(on?0.85:0.1)+'" stroke="'+((bp||hot)?TH.ink:'#fff')+'" stroke-width="'+(bp?1.8:0.6)+
        '" stroke-opacity="'+(on?1:0.12)+'" style="cursor:pointer"><title>'+esc(t.s)+' · '+esc(t.contig||'?')+' from '+esc(t.donor||'?')+
        ' · '+esc(d.lab)+'\nAF in '+(t.af_in==null?'n/a':t.af_in.toFixed(2))+' · AF outside '+(t.af_out==null?'n/a (no site outside the tract)':t.af_out.toFixed(2))+
        ' · '+(t.bp_reads||0)+' breakpoint read(s)'+(t.bf==null?'':'\nlog10 Bayes factor '+t.bf.toFixed(1))+
        '\nclick to highlight this sample</title></circle>';});
    svg+='<text x="'+(padL+plotW/2).toFixed(1)+'" y="'+(H-3)+'" text-anchor="middle" font-size="10" fill="'+TH.mut+'">donor allele fraction OUTSIDE the tract</text>'+
      '<text x="12" y="'+(padT+plotH/2).toFixed(1)+'" text-anchor="middle" font-size="10" fill="'+TH.mut+'" transform="rotate(-90 12 '+(padT+plotH/2).toFixed(1)+')">donor AF inside the tract</text></svg>';
    var legend='<div class="gcv-legend"><span><i class="gcv-ring"></i>a read crosses a breakpoint in cis</span>'+
      '<span class="krk-mut">top-left of the map is what a fixed conversion looks like: donor alleles inside the tract, none outside. Points drift right as donor alleles appear outside it as well, which is what reads arriving from the donor produce. The <b>n/a</b> lane holds tracts with no diagnostic site outside them at all. These two fractions are descriptive: they are what you can go and check in the BAM, while the verdict comes from the model, which reads the bases and their qualities molecule by molecule. Click a point to highlight that sample everywhere.</span></div>';
    // ---- per-tract table: the same evidence as numbers, sortable and searchable ----
    var ord=sel.slice().sort(function(a,b){
      var ka=keyof(a,gconvState.sk),kb=keyof(b,gconvState.sk),d;
      if(typeof ka=='string'||typeof kb=='string'){ka=''+ka;kb=''+kb;d=ka<kb?-1:(ka>kb?1:0);} else d=ka-kb;
      if(!d)d=(b.bp_reads||0)-(a.bp_reads||0)||(a.s<b.s?-1:(a.s>b.s?1:0));
      return gconvState.asc?d:-d;});
    function th(k,lbl,cls,tip){var on=(gconvState.sk===k),ar=on?(gconvState.asc?' &#9650;':' &#9660;'):'';
      return '<th class="gcv-sortable'+(cls?' '+cls:'')+'" data-sk="'+k+'"'+(tip?' title="'+esc(tip)+'"':'')+
        (on?' style="color:var(--accent)"':'')+'>'+lbl+ar+'</th>';}
    var body=ord.slice(0,GCONV_MAXROWS).map(function(t){
      var d=gconvDef(t.verdict), bp=(t.bp_reads||0);
      return '<tr class="gcv-row" data-s="'+esc(t.s)+'"'+(st.hi==t.s?' style="background:'+TH.hl+'"':'')+'>'+
        '<td>'+esc(t.s)+'</td>'+
        '<td>'+gconvLocus(t)+(t.n_don>1?' <span class="gcv-mut" title="'+
           (t.don_call==='resolved'?'the reads pick this relative out of '+t.n_don+' candidates':
            'compatible with '+t.n_don+' relatives; short reads do not carry what would settle it')+
           '">'+(t.don_call==='resolved'?'1 of ':'? of ')+t.n_don+'</span>':'')+'</td>'+
        '<td><span class="gcv-badge" style="background:'+d.c+'" title="'+esc(t.reason||d.tip)+'">'+d.lab+'</span></td>'+
        '<td class="gcv-num"><b'+(t.bf!=null&&t.bf>=3?' class="gcv-bp"':'')+'>'+gconvNum(t.bf,1)+'</b>'+
          (t.bf_null==null?'':' <span class="gcv-mut" title="log10 Bayes factor against no conversion at all. This is where depth, base quality and read linkage show up; the headline number is also limited by how implausible independent substitution is">/ '+t.bf_null.toFixed(0)+' vs none</span>')+'</td>'+
        '<td class="gcv-num">'+(t.n_ev==null?'<span class="gcv-na">n/a</span>':
           (t.n_ev+(t.ev_frac==null?'':' <span class="gcv-mut">'+Math.round(t.ev_frac*100)+'%</span>')))+'</td>'+
        '<td class="gcv-num">'+gconvNum(t.tract_af,2)+'</td>'+
        '<td class="gcv-num">'+gconvNum(t.mismap,3)+'</td>'+
        '<td class="gcv-num">'+(t.start==null?'':fmtpos(t.start))+(t.span==null?'':' <span class="gcv-mut">'+t.span.toLocaleString('en-US')+' bp</span>')+'</td>'+
        '<td class="gcv-num">'+(t.n_sites==null?'':t.n_sites)+(t.n_out==null?'':' <span class="gcv-mut">/ '+t.n_out+' out</span>')+'</td>'+
        '<td class="gcv-num">'+gconvAF(t.af_in,d.c,'no informative depth inside the tract')+'</td>'+
        '<td class="gcv-num">'+gconvAF(t.af_out,'#e0544f','no diagnostic site outside the tract, so boundedness cannot be tested')+'</td>'+
        '<td class="gcv-num"><b'+(bp>0?' class="gcv-bp"':'')+'>'+bp+'</b></td>'+
        '<td class="gcv-num">'+(t.cis_reads==null?'':t.cis_reads)+'</td>'+
        '<td class="gcv-num">'+gconvNum(t.depth,0)+'</td>'+
        '<td><span class="gcv-reason">'+esc(t.reason||'')+'</span></td></tr>';}).join('');
    if(!ord.length)body='<tr><td colspan="15" class="c" style="padding:18px;text-align:center">'+
      (rows.length?'no tract matches the filter.':'no tract in the samples currently in view.')+'</td></tr>';
    host.innerHTML='<div class="gcv-verdicts">'+chips+'</div>'+
      '<div class="gcv-plotscroll">'+svg+'</div>'+legend+
      '<div class="dr-controls">'+
      (hasRep?'<label class="dyn-lab" title="A gene family reports one event once per relationship. Off, each of those rows is listed separately."><input type="checkbox" id="gcvrep"'+(gconvState.oneper?' checked':'')+'> one row per event</label>':'')+
      '<button class="dyn-btn" id="gcvdl" title="Download the tracts listed below as a TSV">'+icon('download')+'download tracts (TSV)</button>'+
      '<span class="dyn-count">'+ord.length+' tract(s)'+(ord.length>GCONV_MAXROWS?' · showing first '+GCONV_MAXROWS+' (download for all)':'')+'</span></div>'+
      '<div class="epitbl-wrap gcv-tablewrap"><table class="epitbl gcv-table"><thead><tr>'+
        th('s','Sample')+th('locus','Locus','','the acceptor locus and the donor its alleles came from')+th('verdict','Verdict')+
        th('bf','BF','gcv-num','log10 Bayes factor for a conversion tract over the best alternative: an independent substitution at the same sites, or reads that arrived from the donor. 3 is decisive')+
        th('n_ev','Samples','gcv-num','how many samples of the cohort carry this event, and what fraction that is. One or two is a finding; nearly all of them means the reference or the aligner, not the isolates')+
        th('tract_af','Carried by','gcv-num','fraction of the reads that carry the tract. Below 1 means either a mixed infection or a third copy of the family contributing unconverted reads; nothing in short reads tells those apart')+
        th('mismap','Donor reads','gcv-num','fraction of reads at this locus the model had to assume came from the donor')+
        th('start','Tract','gcv-num','start position and length of the tract')+
        th('n_sites','Sites','gcv-num','diagnostic sites inside the tract / outside it. No site outside means boundedness cannot be tested')+
        th('af_in','AF in','gcv-num','donor allele fraction inside the tract: near 1 in a clonal sample, intermediate when reads are mismapping')+
        th('af_out','AF out','gcv-num','donor allele fraction at the diagnostic sites OUTSIDE the tract: near 0 when the tract is bounded')+
        th('bp_reads','Bp reads','gcv-num','reads carrying donor alleles one side of a breakpoint and acceptor alleles the other, in cis. A mismapping cannot fake this')+
        th('cis_reads','Cis','gcv-num','reads carrying the donor allele at two or more sites inside the tract')+
        th('depth','Depth','gcv-num','lowest informative depth at any site of the tract')+
        '<th title="which test the verdict rests on">Reason</th></tr></thead><tbody>'+body+'</tbody></table></div>';
    Array.prototype.forEach.call(host.querySelectorAll('.gcv-vchip'),function(c){c.onclick=function(){
      var v=c.getAttribute('data-v'); gconvState.v=(gconvState.v===v)?'':v; draw();};});
    Array.prototype.forEach.call(host.querySelectorAll('.gcv-sortable'),function(h){h.onclick=function(){
      var k=h.getAttribute('data-sk');
      if(gconvState.sk===k)gconvState.asc=!gconvState.asc; else{gconvState.sk=k; gconvState.asc=(k=='s'||k=='locus'||k=='verdict');}
      draw();};});
    Array.prototype.forEach.call(host.querySelectorAll('[data-s]'),function(e){e.onclick=function(){setHi(e.getAttribute('data-s'));};});
    var rb=el('gcvrep'); if(rb)rb.onchange=function(){gconvState.oneper=this.checked; draw();};
    var db=el('gcvdl'); if(db)db.onclick=function(){
      // Laid out exactly like the cohort TSV on disk: the tool's own columns, then the ones the
      // cohort pass appends. The verdict column is the sample's own and cohort_verdict is what
      // the rest of the cohort made of it, which is the pair the file itself carries.
      var hdr=['sample','pair_id','contig','donor','verdict','reason','start','end','span_bp',
               'don_start','don_end',
               'post_conv','log10_bf','log10_bf_vs_null','tract_af','mismap_frac','mut_rate','start_ci','end_ci',
               'n_sites','n_sites_outside','n_undetermined','donor_af_in','donor_af_outside',
               'min_depth','cis_reads','breakpoint_reads','donor_only_reads',
               'n_derived','n_ancestral','n_unpolarised','donor_swap_af','donor_swap_af_outside',
               'locus_cn','expected_af',
               'genes','n_syn','n_nonsyn','aa_changes',
               'event_id','event_samples','event_frac','cohort_verdict','cohort_mismap','cohort_bf_median',
               'donor_rank','n_donors','donor_margin','donor_call','is_representative'];
      var lines=[hdr.join('\t')];
      ord.forEach(function(t){lines.push([t.s,t.pair,t.contig,t.donor,(t.sample_verdict||t.verdict),t.reason,t.start,t.end,t.span,
        t.don_start,t.don_end,
        t.post,t.bf,t.bf_null,t.tract_af,t.mismap,t.mut_rate,t.start_ci,t.end_ci,
        t.n_sites,t.n_out,t.n_undet,t.af_in,t.af_out,t.depth,t.cis_reads,t.bp_reads,t.donor_only,
        t.n_derived,t.n_ancestral,t.n_unpol,t.donor_swap,t.donor_swap_out,
        t.locus_cn,t.expected_af,
        t.genes,t.n_syn,t.n_nonsyn,t.aa_changes,
        t.event,t.n_ev,t.ev_frac,t.verdict,t.co_mismap,t.co_bf,
        t.don_rank,t.n_don,t.don_margin,t.don_call,t.rep].map(function(x){return x==null?'':x;}).join('\t'));});
      dl(lines.join('\n')+'\n','gene_conversion.tsv','text/tab-separated-values');};
    if(cap){var nsamp={},ncall=0,nev={}; rows.forEach(function(t){nsamp[t.s]=1;
        if(t.event!=null&&t.event!=='')nev[t.event]=1;
        if(t.verdict==='gene_conversion'&&(t.rep==null||t.rep===1))ncall++;});
      var evtxt=Object.keys(nev).length?(Object.keys(nev).length+' event(s) over '+rows.length+' row(s)')
                                       :(rows.length+' candidate tract(s)');
      cap.innerHTML=evtxt+' in '+Object.keys(nsamp).length+' sample(s) of the cohort in view &#183; <b'+
        (ncall?' class="sc-sig"':'')+'>'+ncall+' called gene conversion</b>, '+(counts.ambiguous||0)+' ambiguous, '+(counts.mismapping||0)+
        ' mismapping, '+(counts.coverage_shift||0)+' coverage shift, '+(counts.reference_artifact||0)+' reference artifact, '+(counts.reference_derived||0)+' reference derived, '+(counts.reciprocal_exchange||0)+' reciprocal exchange &#183; <b>'+nbp+'</b> supported by a read crossing a breakpoint in cis. <span class="krk-mut">These loci are repeats, so they are excluded from variant calling and the consensus by design: a tract will not appear in the SNP matrix, and that is expected. Breakpoints are located to diagnostic-site resolution, not to the base, and each paralog pair is judged independently, so one tract can be reported against more than one donor. Candidates to inspect, not confirmed events.</span>'+gconvSettings();}
  }
  draw();
}

