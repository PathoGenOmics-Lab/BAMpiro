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

