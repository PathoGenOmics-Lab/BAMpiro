function drGColor(gn){ return (gn===1||gn===2)?'#dc2626':(gn===3?'#d97706':((gn===4||gn===5)?'#94a3b8':'#b8c2cf')); }
function drGInk(gn){ return (gn===1||gn===2||gn===3)?'#fff':'#1f2a37'; }   // dark ink on the pale grey (4-5 / unknown) badges so the label stays legible
// The last branch used to be reached two ways: no mutation at all, and a mutation whose grade is
// not a WHO number. Both rendered as an empty cell, so a detected variant the catalogue does not
// grade looked exactly like a clean drug. Catalogue v1.0.0 had 6,056 such rows and ten of them are
// WHO grade 1-2 in v1.0.2. An ungraded call now gets a mark of its own; only a genuine no-call is blank.
function drStatus(gns){ var r=false,u=false,x=false,n=false; for(var i=0;i<gns.length;i++){var g=gns[i]; if(g===1||g===2)r=true; else if(g===3)u=true; else if(g===4||g===5)n=true; else x=true;} return r?{t:'R',c:'#dc2626'}:(u?{t:'?',c:'#d97706'}:(x?{t:'!',c:'#7c6f9f'}:(n?{t:'&#183;',c:'#94a3b8'}:{t:'',c:'var(--track)'}))); }
// ---- Kraken2: one row per sample, worst first; the flagged ones by default ----
var krkView={all:false};
var KRK_TAXPAL=['#e0544f','#e0a11f','#8a63c9','#d06fae','#26a0a0','#c98a3b','#3b7dd8','#7a8794'];
var KRK_STATE={off:['other organism','bad'],mixed:['mixed','warn'],uncl:['many unclassified','warn'],ok:['clean','good']};
function renderKraken(){
  var host=el('kraken_body'), sec=el('kraken'); if(!host)return;
  var K=R.kraken;
  if(!(K&&K.samples&&K.samples.length)){ if(sec)sec.style.display='none'; var nv=el('nav-kraken'); if(nv)nv.style.display='none'; return; }
  var tgt=krkTarget();
  var all=K.samples.slice().sort(function(a,b){return (krkPurity(a)-krkPurity(b))||(b.unclassified-a.unclassified)||(a.s<b.s?-1:1);});
  var flagged=all.filter(function(k){return krkState(k)!=='ok';});
  var showAll=krkView.all||!flagged.length, rows=showAll?all:flagged;
  // one colour per taxon across every row: the target in the pass colour, then the most abundant others
  var agg={}; rows.forEach(function(k){(k.top||[]).forEach(function(t){if(t.name!==tgt)agg[t.name]=(agg[t.name]||0)+t.pct;});});
  var others=Object.keys(agg).sort(function(a,b){return agg[b]-agg[a];}).slice(0,KRK_TAXPAL.length);
  var tcol={}; tcol[tgt]='var(--pass)'; others.forEach(function(t,i){tcol[t]=KRK_TAXPAL[i];});
  function comp(k){var used=0,h='';
    (k.top||[]).forEach(function(t){var c=tcol[t.name]; if(!c||!(t.pct>0))return; used+=t.pct;
      h+='<span style="width:'+t.pct.toFixed(2)+'%;background:'+c+'" title="'+esc(t.name)+' '+t.pct.toFixed(1)+'% of reads"></span>';});
    var other=Math.max(0,(k.classified||0)-used); if(other>0.05)h+='<span style="width:'+other.toFixed(2)+'%;background:#9aa7b6" title="other classified '+other.toFixed(1)+'%"></span>';
    if(k.unclassified>0.05)h+='<span style="width:'+k.unclassified.toFixed(2)+'%;background:var(--track)" title="unclassified '+k.unclassified.toFixed(1)+'%"></span>';
    return '<div class="krk-bar">'+h+'</div>';}
  function taxon(t){return t?'<i>'+esc(t.name)+'</i> <span class="krk-mut">'+t.pct.toFixed(1)+'%</span>':'<span class="krk-mut">-</span>';}
  var body=rows.map(function(k){var st0=krkState(k),lab=KRK_STATE[st0],p=krkPurity(k);
    var next=(k.secondary&&k.secondary.pct>=0.5)?k.secondary:null;
    return '<tr class="hit" data-s="'+esc(k.s)+'"'+(st.hi==k.s?' style="background:'+TH.hl+'"':'')+'><td class="s">'+esc(k.s)+'</td>'+
      '<td style="text-align:left"><span class="tag '+lab[1]+'">'+lab[0]+'</span></td>'+
      '<td class="num"><span class="krk-pur"><span class="krk-purbar"><span style="width:'+Math.max(0,Math.min(100,p)).toFixed(1)+'%;background:'+(st0==='off'?'var(--fail)':(st0==='mixed'?'var(--warn)':'var(--pass)'))+'"></span></span>'+p.toFixed(1)+'</span></td>'+
      '<td style="text-align:left">'+taxon(k.primary)+'</td>'+
      '<td style="text-align:left">'+taxon(next)+'</td>'+
      '<td'+(k.unclassified>KRK_UNCL?' style="color:var(--warn);font-weight:600"':'')+'>'+k.unclassified.toFixed(1)+'</td>'+
      '<td>'+(k.runs||1)+'</td>'+
      '<td class="krk-barcell">'+comp(k)+'</td></tr>';}).join('');
  var legend='<div class="krk-legend">'+[tgt].concat(others).map(function(t){return '<span><i style="background:'+tcol[t]+'"></i>'+esc(t)+'</span>';}).join('')+
    '<span><i style="background:#9aa7b6"></i>other classified</span><span><i style="background:var(--track)"></i>unclassified</span></div>';
  host.innerHTML='<div class="dr-controls">'+
      '<span class="seg" id="krkseg"><button data-v="f"'+(!showAll?' class="on"':'')+(flagged.length?'':' disabled')+'>flagged <b>'+flagged.length+'</b></button><button data-v="a"'+(showAll?' class="on"':'')+'>all <b>'+all.length+'</b></button></span>'+
      '<span class="c">sorted by the share of classified reads in <i>'+esc(tgt)+'</i>, lowest first; below '+KRK_PURE+'% is mixed, below '+KRK_OFF+'% another organism</span></div>'+
    legend+
    '<div class="gtable krk-table"><table class="krktable"><thead><tr><th class="s">Sample</th><th style="text-align:left">Status</th><th title="share of the classified reads in '+esc(tgt)+', the cohort target">In target %</th><th style="text-align:left">Dominant taxon</th><th style="text-align:left">Next taxon</th><th>Unclass. %</th><th title="Kraken2 reports merged for this sample">Runs</th><th style="text-align:left">Composition of all reads</th></tr></thead><tbody>'+body+'</tbody></table></div>';
  Array.prototype.forEach.call(host.querySelectorAll('#krkseg button'),function(b){b.onclick=function(){krkView.all=(b.getAttribute('data-v')==='a');renderKraken();};});
  Array.prototype.forEach.call(host.querySelectorAll('tr[data-s]'),function(e){e.onclick=function(){setHi(e.getAttribute('data-s'));};});
}
// ---- Drug resistance: mutations first (grouped across samples), then the sample x drug matrix and every call ----
var drState={q:'',view:'mut',g:{r:true,u:false,n:false,x:false},open:{}};
var DR_GRADES=[['r','1–2 associated','#dc2626'],['u','3 uncertain','#d97706'],['n','4–5 not associated','#94a3b8'],['x','ungraded','#7c6f9f']];
function drGClass(gn){return (gn===1||gn===2)?'r':(gn===3?'u':((gn===4||gn===5)?'n':'x'));}
function renderDrug(){
  var host=el('drug_body'), sec=el('drug'); if(!host)return;
  var D=R.dr;
  if(!(D&&D.calls&&D.calls.length)){ if(sec)sec.style.display='none'; var nv=el('nav-drug'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  var sum=drSummary(), calls=D.calls, drugs=D.drugs;
  var linOf={}; R.samples.forEach(function(s){linOf[s.s]=untyped(s.lineage)?null:linMain(s.lineage);});
  var wideKey={}; sum.wide.forEach(function(m){wideKey[m.key]=m.wide;});
  function gradeOn(c){return !!drState.g[drGClass(c.gn)];}
  function qOk(c){var q=drState.q.toLowerCase(); return !q||(c.s.toLowerCase().indexOf(q)>=0)||(c.drug.toLowerCase().indexOf(q)>=0)||(c.gene.toLowerCase().indexOf(q)>=0)||(c.mutation.toLowerCase().indexOf(q)>=0);}
  var gradeBtns=DR_GRADES.map(function(g){var n=calls.filter(function(c){return drGClass(c.gn)===g[0];}).length;
    return '<button data-g="'+g[0]+'"'+(drState.g[g[0]]?' class="on"':'')+(n?'':' disabled')+'><i style="background:'+g[2]+'"></i>'+g[1]+' <b>'+n+'</b></button>';}).join('');
  host.innerHTML='<div class="dr-controls">'+
      '<span class="seg" id="drview"><button data-v="mut"'+(drState.view==='mut'?' class="on"':'')+'>mutations</button><button data-v="mx"'+(drState.view==='mx'?' class="on"':'')+'>samples &#215; drugs</button><button data-v="calls"'+(drState.view==='calls'?' class="on"':'')+'>every call</button></span>'+
      '<span class="seg dr-grades" id="drgrades" title="WHO confidence grades to show">'+gradeBtns+'</span>'+
      '<input id="drq" class="dyn-search" type="search" placeholder="filter by sample / drug / gene / mutation..." value="'+esc(drState.q)+'">'+
      '<button class="dyn-btn" id="drdl" title="Download every resistance call as a TSV">'+icon('download')+'download calls (TSV)</button>'+
      '<span class="dyn-count" id="drcount"></span></div><div id="drview_body"></div>';
  function mutView(){
    var agg={},order=[];
    calls.forEach(function(c){if(!gradeOn(c)||!qOk(c))return;var k=c.drug+'|'+c.gene+'|'+c.mutation,m=agg[k];
      if(!m){m=agg[k]={key:k,drug:c.drug,gene:c.gene,mutation:c.mutation,gn:c.gn,grade:c.grade,samples:[],afs:[]};order.push(k);}
      if(m.samples.indexOf(c.s)<0)m.samples.push(c.s); if(c.af!=null)m.afs.push(c.af);});   // a sample merged from several runs can list a mutation once per run
    var list=order.map(function(k){return agg[k];});
    list.sort(function(a,b){return ((wideKey[a.key]?1:0)-(wideKey[b.key]?1:0))||((a.gn||9)-(b.gn||9))||(b.samples.length-a.samples.length)||(a.key<b.key?-1:1);});
    el('drcount').textContent=list.length+' mutation(s) in '+calls.filter(function(c){return gradeOn(c)&&qOk(c);}).length+' call(s)';
    if(!list.length)return '<div class="nd pad">No mutation at the selected grades'+(drState.q?' matches the filter':'')+'.</div>';
    var nLin=sum.nLin;
    return '<div class="epitbl-wrap"><table class="epitbl drmut"><thead><tr><th>Drug</th><th>Gene</th><th>Mutation (H37Rv)</th><th>WHO grade</th><th class="num">Samples</th><th>By lineage</th><th class="num">AF (median)</th><th></th></tr></thead><tbody>'+
      list.map(function(m){var per={}; m.samples.forEach(function(s){var l=linOf[s]||'untyped';per[l]=(per[l]||0)+1;});
        var wl=wideKey[m.key], afm=m.afs.length?_median(m.afs):null, open=!!drState.open[m.key];
        var lin=Object.keys(per).sort(function(a,b){return per[b]-per[a];}).map(function(l){return esc(l)+' '+per[l]+(nLin[l]?'/'+nLin[l]:'');}).join(' &#183; ');
        return '<tr class="dr-mrow'+(open?' open':'')+'" data-k="'+esc(m.key)+'"><td><b>'+esc(m.drug)+'</b></td><td>'+esc(m.gene)+geneRvTag(m.gene)+'</td><td class="epitbl-r">'+esc(m.mutation)+'</td>'+
          '<td><span class="dr-badge" style="background:'+drGColor(m.gn)+';color:'+drGInk(m.gn)+'">'+esc(m.grade||'?')+'</span></td>'+
          '<td class="num"><b>'+m.samples.length+'</b></td><td class="dr-lin">'+lin+'</td><td class="num">'+(afm==null?'':afm.toFixed(2))+'</td>'+
          '<td>'+(wl?'<span class="tag neu" title="carried by at least '+Math.round(DR_WIDE*100)+'% of the '+esc(wl)+' samples: the lineage\'s own, not something that arose here">'+esc(wl)+' marker</span>':'')+'</td></tr>'+
          (open?'<tr class="dr-carriers"><td colspan="8">'+m.samples.slice().sort().map(function(s){return '<span class="dr-car" data-s="'+esc(s)+'">'+esc(s)+'</span>';}).join('')+'</td></tr>':'');}).join('')+
      '</tbody></table></div><div class="krk-mut" style="padding:6px 2px">Click a mutation to list the samples that carry it. A lineage marker is a mutation at least '+Math.round(DR_WIDE*100)+'% of a lineage carries.</div>';
  }
  function mxView(){
    var cell={}, cmut={}, dset={};
    calls.forEach(function(c){ if(!gradeOn(c)||!qOk(c))return; var ds=(c.dr&&c.dr.length)?c.dr:[c.drug]; for(var i=0;i<ds.length;i++){ var d=ds[i]; dset[d]=1;
      (cell[c.s]=cell[c.s]||{}); (cell[c.s][d]=cell[c.s][d]||[]).push(c.gn); (cmut[c.s]=cmut[c.s]||{}); (cmut[c.s][d]=cmut[c.s][d]||[]).push(c); } });
    var cols=drugs.filter(function(d){return dset[d];});
    // samples with a mutation their lineage does not share first, then any call, then the rest
    var rank=function(s){return sum.carriers[s]?0:(cell[s]?1:2);};
    var rows=(D.samples||[]).filter(function(s){return cell[s];}).sort(function(a,b){return rank(a)-rank(b)||(a<b?-1:1);});
    el('drcount').textContent=rows.length+' sample(s) with a call at the selected grades';
    if(!rows.length||!cols.length)return '<div class="nd pad">No call at the selected grades.</div>';
    return '<div class="dr-legend"><span><i style="background:#dc2626"></i>R: grade 1&#8211;2</span><span><i style="background:#d97706"></i>?: grade 3</span><span><i style="background:#94a3b8"></i>&#183;: grade 4&#8211;5</span><span><i style="background:#7c6f9f"></i>!: ungraded</span>'+
      '<span class="c">worst grade per drug; blank = no call. Samples with a mutation their lineage does not share come first.</span></div>'+
      '<div class="dr-mxwrap"><table class="drmx"><thead><tr><th class="dr-corner">sample \\ drug</th>'+
      cols.map(function(dr){return '<th class="dr-hcell" title="'+esc(dr)+'"><span class="dr-h">'+esc(dr)+'</span></th>';}).join('')+'</tr></thead><tbody>'+
      rows.map(function(s){ return '<tr><th class="dr-row'+(sum.carriers[s]?' acq':'')+'" title="'+esc(s)+(linOf[s]?' · '+esc(linOf[s]):'')+'">'+esc(s)+'</th>'+cols.map(function(dr){
        var gns=(cell[s]||{})[dr];
        if(!gns) return '<td class="drmx-cell" title="'+esc(s)+' &#183; '+esc(dr)+': no call"></td>';
        var st0=drStatus(gns);
        var muts=((cmut[s]||{})[dr]||[]).map(function(c){return c.gene+' '+c.mutation+(c.gn?(' (WHO '+c.gn+')'):'');}).join('; ');
        return '<td class="drmx-cell" style="background:'+st0.c+'" title="'+esc(s)+' &#183; '+esc(dr)+' &#8212; '+esc(muts)+'"><b>'+st0.t+'</b></td>';
      }).join('')+'</tr>'; }).join('')+'</tbody></table></div>';
  }
  function callView(){
    var rows=calls.filter(function(c){return gradeOn(c)&&qOk(c);});
    rows=rows.slice().sort(function(a,b){ if(a.s!==b.s) return a.s<b.s?-1:1; return (a.gn||9)-(b.gn||9); });
    el('drcount').textContent=rows.length+' call(s)'+(rows.length>600?' · showing the first 600 (download for all)':'');
    var h='<div class="epitbl-wrap"><table class="epitbl"><thead><tr><th>Sample</th><th>Drug</th><th>Gene</th><th>Mutation (H37Rv)</th><th>WHO grade</th><th>AF</th><th>DP</th></tr></thead><tbody>';
    if(!rows.length) h+='<tr><td colspan="7" class="c" style="padding:18px;text-align:center">no call matches the filter.</td></tr>';
    h+=rows.slice(0,600).map(function(c){
      return '<tr><td>'+esc(c.s)+'</td><td><b>'+esc(c.drug)+'</b></td><td>'+esc(c.gene)+geneRvTag(c.gene)+'</td><td class="epitbl-r">'+esc(c.mutation)+'</td>'+
        '<td><span class="dr-badge" style="background:'+drGColor(c.gn)+';color:'+drGInk(c.gn)+'" title="'+esc(c.marker||'')+'">'+esc(c.grade||'?')+'</span></td>'+
        '<td>'+(c.af==null?'':c.af.toFixed(2))+'</td><td>'+(c.dp==null?'':c.dp)+'</td></tr>';
    }).join('')+'</tbody></table></div>';
    return h;
  }
  function draw(){
    var vb=el('drview_body'); vb.innerHTML=drState.view==='mx'?mxView():(drState.view==='calls'?callView():mutView());
    Array.prototype.forEach.call(vb.querySelectorAll('.dr-mrow'),function(tr){tr.onclick=function(){var k=tr.getAttribute('data-k');drState.open[k]=!drState.open[k];draw();};});
    Array.prototype.forEach.call(vb.querySelectorAll('.dr-car'),function(sp){sp.onclick=function(e){e.stopPropagation();openDetail(sp.getAttribute('data-s'));};});
  }
  el('drq').oninput=function(){ drState.q=this.value; draw(); };
  Array.prototype.forEach.call(host.querySelectorAll('#drview button'),function(b){b.onclick=function(){drState.view=b.getAttribute('data-v');
    Array.prototype.forEach.call(host.querySelectorAll('#drview button'),function(x){x.classList.toggle('on',x===b);});draw();};});
  Array.prototype.forEach.call(host.querySelectorAll('#drgrades button'),function(b){b.onclick=function(){var g=b.getAttribute('data-g');drState.g[g]=!drState.g[g];b.classList.toggle('on',drState.g[g]);draw();};});
  el('drdl').onclick=function(){
    var hdr=['sample','drug','gene','mutation_h37rv','who_grade','marker','af','dp','lineage_marker'];
    var lines=[hdr.join('\t')];
    calls.forEach(function(c){ lines.push([c.s,c.drug,c.gene,c.mutation,c.grade,c.marker,(c.af==null?'':c.af),(c.dp==null?'':c.dp),(wideKey[c.drug+'|'+c.gene+'|'+c.mutation]||'')].join('\t')); });
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
var gconvState={q:'',v:'',sk:'verdict',asc:true,oneper:true};   // called events first, and within a verdict the strongest (BF) first
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
              tip:'present in nearly every sample mapped to the same reference. The reference being wrong here, or the aligner doing this to everybody, explains that more simply than the same conversion arising in every isolate. In a CLONAL cohort it may instead be shared ancestry, which recurrence alone cannot distinguish. Only a cohort can make this call at all'},
             {k:'reciprocal_exchange',lab:'reciprocal exchange',c:'#c77d3a',r:6,
              tip:'the donor carries the ACCEPTOR\'s bases over the same stretch, so both copies changed. That is an exchange between them rather than one being overwritten, and gene conversion is non-reciprocal by definition'},
             {k:'reference_derived',lab:'reference derived',c:'#4a90b8',r:5,
              tip:'an outgroup says the REFERENCE carries the derived base over this stretch and the reads carry the ancestral one. The sample changed nothing; the finding belongs to the reference. Without an outgroup this is the same picture as a conversion'},
             // The cohort step sets this for every tract of a sample that calls tracts at too many loci
             // for them to be local events. It was missing here, so the most common verdict of a run
             // had no chip and every one of its rows read "verdict not recognised".
             {k:'divergent_sample',lab:'divergent sample',c:'#6b7280',r:7,
              tip:'the sample calls tracts at too many loci for them to be local events: its genome differs from the reference as a whole (often a sample mapped to the wrong lineage), so none of its tracts is read as conversion'}];
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
      if(d)return gconvState.asc?d:-d;
      return ((b.bf==null?-1e9:b.bf)-(a.bf==null?-1e9:a.bf))||((b.bp_reads||0)-(a.bp_reads||0))||(a.s<b.s?-1:(a.s>b.s?1:0));});
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
        '<td><span class="gcv-reason" title="'+esc(t.reason||'')+'">'+esc(t.reason||'')+'</span></td></tr>';}).join('');
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
        th('n_ev','Samples','gcv-num','how many samples carry this event, and what fraction that is of the samples mapped to the same reference. One or two is a finding; nearly all of them means the reference or the aligner, not the isolates')+
        th('tract_af','Carried by','gcv-num','fraction of the reads that carry the tract. Below 1 means either a mixed infection or a third copy of the family contributing unconverted reads; nothing in short reads tells those apart. The model does not go below 20%, so a thinner share is shown as 20%: AF in is what the reads themselves carry')+
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

