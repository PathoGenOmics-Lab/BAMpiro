var EPICOL={A:'#2f6fed',B:'#e6893a'};
var epiState={dir:'all',q:'',minr:null,conf:'all',view:'cards',tsort:{k:'q',asc:true},showAll:false};
var epiZoom=250;   // epistasis pair-card width in px (zoom slider), mirrors the dynamics cards
function epiMiniChart(p){
  var times=p.times||[], A=p.trajA||[], B=p.trajB||[], n=times.length;
  var W=250,H=150,ml=30,mr=12,mt=10,mb=26,pw=W-ml-mr,ph=H-mt-mb;
  function X(i){ return ml+(n<=1?pw/2:(i/(n-1))*pw); }
  function Y(a){ return mt+(1-a)*ph; }
  var svg='<svg viewBox="0 0 '+W+' '+H+'" width="100%" style="display:block"><title>Two allele-frequency trajectories over time; parallel lines = concordant, mirrored = discordant. Hover a point for its value.</title>';
  [0,0.5,1].forEach(function(a){ svg+='<line x1="'+ml+'" y1="'+Y(a).toFixed(1)+'" x2="'+(W-mr)+'" y2="'+Y(a).toFixed(1)+'" stroke="'+TH.grid+'" stroke-width="1"/><text x="'+(ml-6)+'" y="'+(Y(a)+3.5).toFixed(1)+'" text-anchor="end" font-size="10" fill="'+TH.mut+'">'+a.toFixed(1)+'</text>'; });
  svg+='<line x1="'+ml+'" y1="'+mt+'" x2="'+ml+'" y2="'+(mt+ph).toFixed(1)+'" stroke="'+TH.axis+'" stroke-width="1"/>';   // left value axis + ticks (matches the dynamics cards)
  [0,0.5,1].forEach(function(a){ svg+='<line x1="'+(ml-3)+'" y1="'+Y(a).toFixed(1)+'" x2="'+ml+'" y2="'+Y(a).toFixed(1)+'" stroke="'+TH.axis+'" stroke-width="1"/>'; });
  [[A,EPICOL.A],[B,EPICOL.B]].forEach(function(pr){ var t=pr[0],c=pr[1],lastI=t.length-1;
    var pts=t.map(function(a,i){return X(i).toFixed(1)+','+Y(a).toFixed(1);}).join(' ');
    svg+='<polyline points="'+pts+'" fill="none" stroke="'+c+'" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/>';
    t.forEach(function(a,i){ if(i===lastI)return; var cx=X(i).toFixed(1), cy=Y(a).toFixed(1); svg+='<circle cx="'+cx+'" cy="'+cy+'" r="3" fill="'+TH.panel+'"/><circle cx="'+cx+'" cy="'+cy+'" r="2.1" fill="'+c+'"/>'; });   // panel-haloed beads
    var ex=X(lastI).toFixed(1), ey=Y(t[lastI]).toFixed(1);
    svg+='<circle cx="'+ex+'" cy="'+ey+'" r="4" fill="'+TH.panel+'" stroke="'+c+'" stroke-width="2.2"/><circle cx="'+ex+'" cy="'+ey+'" r="1.4" fill="'+c+'"/>';   // hollow last-value ring
    t.forEach(function(a,i){ svg+='<circle cx="'+X(i).toFixed(1)+'" cy="'+Y(a).toFixed(1)+'" r="6" fill="transparent"><title>t='+esc(times[i]==null?i:times[i])+'  AF='+a.toFixed(3)+'</title></circle>'; });   // invisible hit targets keep tooltips
  });
  times.forEach(function(t,i){ svg+='<text x="'+X(i).toFixed(1)+'" y="'+(H-8)+'" text-anchor="middle" font-size="10.5" fill="'+TH.mut+'">'+esc(t==null?i:t)+'</text>'; });
  svg+='</svg>';
  return svg;
}
function epiColor(r){
  if(r==null) return '#f3f5f8';
  var a=Math.min(1,Math.abs(r));
  if(r>=0) return 'rgba(47,143,91,'+(0.10+a*0.82).toFixed(2)+')';   // concordant (green)
  return 'rgba(162,74,143,'+(0.10+a*0.82).toFixed(2)+')';           // discordant (purple)
}
function epiFmtR(v){ return (v<0?'-':(v>0?'+':''))+Math.abs(v).toFixed(2).replace(/^0/,''); }
function epiNodeLbl(nd){ return (nd.gene||'(intergenic)')+' '+String(nd.pos).split(':').pop(); }
function renderEpistasis(){
  var host=el('epi_body'), sec=el('epistasis'); if(!host)return;
  var E=R.epistasis;
  if(!(E&&E.pairs&&E.pairs.length)){ if(sec)sec.style.display='none'; var nv=el('nav-epi'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  if(epiState.minr==null) epiState.minr=E.min_r||0.8;
  var DIRS=[['all','all pairs'],['concordant','same dynamics'],['discordant','opposite dynamics']];
  var CONF=[['all','all pairs'],['mod','moderate or better (permutation p &#8804; 0.05)'],['strong','strong only (FDR q &#8804; 0.05)']];
  var VIEWS=[['cards','cards'],['matrix','matrix'],['table','table']];
  host.innerHTML=
    '<div class="epi-views">'+VIEWS.map(function(v){return '<button class="epi-viewbtn'+(epiState.view===v[0]?' on':'')+'" data-v="'+v[0]+'">'+v[1]+'</button>';}).join('')+'</div>'+
    '<div class="epi-controls">'+
      '<input id="epiq" class="dyn-search" type="search" title="Filter the pairs by gene name or position" placeholder="filter by gene / position..." value="'+esc(epiState.q)+'">'+
      '<span class="epi-flabel">dynamics</span>'+
      DIRS.map(function(d){return '<button class="dyn-btn epi-dirbtn'+(epiState.dir===d[0]?' on':'')+'" data-d="'+d[0]+'" title="Show '+d[1]+'">'+d[0]+'</button>';}).join('')+
      '<span class="epi-flabel">confidence</span>'+
      CONF.map(function(c){return '<button class="dyn-btn epi-confbtn'+(epiState.conf===c[0]?' on':'')+'" data-c="'+c[0]+'" title="'+c[1]+'">'+c[0]+'</button>';}).join('')+
      '<label class="dyn-zoom" title="Minimum |Pearson r| for a pair to be shown (cards / table)"><span>|r| &#8805;</span><input type="range" id="epir" min="'+(E.min_r||0.8)+'" max="0.99" step="0.01" value="'+epiState.minr+'"><b id="epirv">'+epiState.minr.toFixed(2)+'</b></label>'+
      '<label class="dyn-zoom" title="Resize the pair cards - drag left to fit more per row"><span>'+icon('search','sort')+'&#8211;/+</span><input type="range" id="epizoom" min="165" max="360" step="5" value="'+epiZoom+'"></label>'+
      '<button class="dyn-btn showall-btn" id="epiShowAll" title="Show every reported pair (clear the direction / confidence / |r| filters)">show all</button>'+
      '<span class="dyn-count" id="epicount"></span></div>'+
    '<div class="epi-legend">'+
      '<span title="The two variants rise and fall together across the series - candidate linkage or co-selection."><i style="background:#2f8f5b"></i>concordant / same dynamics <span class="infoi">i</span></span>'+
      '<span title="One variant rises as the other falls across the series - competing lineages / clonal interference."><i style="background:#a24a8f"></i>discordant / opposite dynamics <span class="infoi">i</span></span>'+
      '<span title="Permutation p-value: how often shuffling the timepoints of one trajectory reaches this |r| by chance. q = Benjamini-Hochberg FDR across every reported pair. strong = q &#8804; 0.05, moderate = p &#8804; 0.05, weak otherwise. With few timepoints a single series cannot beat ~1/n! by chance, so a pattern RECURRING across independent series is what drives a pair to strong."><b>strong</b> &#183; moderate &#183; weak = permutation p &amp; FDR q <span class="infoi">i</span></span>'+
      '<span class="c">mean Pearson r within a series (&#8805; '+E.min_points+' timepoints), across '+E.n_series+' series; '+E.perm+' permutations.'+
        ((E.vars_capped||E.pairs_capped)?(' <b>Large cohort:</b> '+(E.vars_capped?('top '+E.max_vars+' most-variable variants/series (of &#8804;'+E.n_vars_max+')'):'')+(E.vars_capped&&E.pairs_capped?', ':'')+(E.pairs_capped?(E.max_pairs_perm+' strongest of '+E.n_candidates+' candidate pairs tested'):'')+'.'):'')+'</span>'+
    '</div>'+
    '<div class="dyn-grid" id="epi-cards"></div>'+
    '<div id="epi-matrix"></div>'+
    '<div id="epi-table"></div>';
  function posn(x){ return String(x).split(':').pop(); }
  function confok(p){ if(epiState.conf==='strong')return p.tier==='strong'; if(epiState.conf==='mod')return p.tier==='strong'||p.tier==='moderate'; return true; }
  function epiFilter(){
    var q=epiState.q.toLowerCase();
    return E.pairs.filter(function(p){
      if(Math.abs(p.r)<epiState.minr) return false;
      if(epiState.dir!=='all'&&p.direction!==epiState.dir) return false;
      if(!confok(p)) return false;
      if(q && !((p.geneA&&p.geneA.toLowerCase().indexOf(q)>=0)||(p.geneB&&p.geneB.toLowerCase().indexOf(q)>=0)||String(p.posA).indexOf(q)>=0||String(p.posB).indexOf(q)>=0)) return false;
      return true;
    });
  }
  function epiCards(list){
    var grid=el('epi-cards');
    grid.style.setProperty('--dyncw', epiZoom+'px');
    if(!list.length){ grid.innerHTML='<div class="dyn-empty" style="grid-column:1/-1">&#128204; no variant pair matches the current filter.</div>'; return; }
    var show=epiState.showAll?list:list.slice(0,12);   // collapsed shows the top pairs; 'show all' switch reveals them all
    grid.innerHTML=show.map(function(p){
      var arrow=p.direction==='concordant'?'&#8596;':'&#8646;';
      var rec=(p.n>1)?('<span class="epi-recur" title="seen in '+p.n+' independent series'+(p.consistent?' with the same sign - recurrent':'')+'">&#8635; '+p.n+' series</span>'):('<span title="from a single series">series '+esc(p.group)+'</span>');
      return '<div class="epi-card '+p.direction+'">'+
        '<div class="epi-card-h"><span class="epi-badge '+p.direction+'">r = '+(p.r>0?'+':'')+p.r.toFixed(2)+'</span><span class="epi-tier epi-'+p.tier+'" title="permutation p = '+p.p+', FDR q = '+p.q+'">'+p.tier+'</span></div>'+
        '<div class="epi-pair"><span style="color:'+EPICOL.A+'"><b>'+esc(p.geneA||'(intergenic)')+'</b>'+geneRvTag(p.geneA)+' '+posn(p.posA)+(p.aaA?(' '+aaDual(p.aaA,p.aaA_h37rv)):'')+'</span><span class="epi-vs">'+arrow+'</span><span style="color:'+EPICOL.B+'"><b>'+esc(p.geneB||'(intergenic)')+'</b>'+geneRvTag(p.geneB)+' '+posn(p.posB)+(p.aaB?(' '+aaDual(p.aaB,p.aaB_h37rv)):'')+'</span></div>'+
        epiMiniChart(p)+
        '<div class="epi-card-f"><span title="permutation p-value / Benjamini-Hochberg FDR q-value">p '+p.p.toFixed(3)+' &#183; q '+p.q.toFixed(3)+'</span>'+rec+'</div>'+
      '</div>';
    }).join('');
  }
  function epiMatrix(){
    var box=el('epi-matrix'), M=E.matrix;
    if(!(M&&M.nodes&&M.nodes.length)){ box.innerHTML='<div class="dyn-empty">not enough correlated variants for a matrix.</div>'; return; }
    var nodes=M.nodes, N=nodes.length, showVal=N<=18;
    function cell(i,j){ if(i===j)return 1; var lo=Math.min(i,j),hi=Math.max(i,j); var v=M.cells[lo+','+hi]; return (v==null)?null:v; }
    var h='<div class="epimx-note">'+N+' variant(s)'+(M.truncated?(' (top '+N+' of '+M.total_nodes+' by connectivity)'):'')+' &#183; every pairwise mean r, including sub-threshold cells &#183; hover a cell for the pair.</div>';
    h+='<div class="epimx-wrap"><table class="epimx"><thead><tr><th class="epimx-corner"></th>';
    nodes.forEach(function(nd){ h+='<th class="epimx-hcell" title="'+esc(epiNodeLbl(nd))+(nd.aa?(' '+esc(nd.aa)):'')+'"><span class="epimx-h">'+esc(epiNodeLbl(nd))+'</span></th>'; });
    h+='</tr></thead><tbody>';
    nodes.forEach(function(nd,i){
      h+='<tr><th class="epimx-row" title="'+esc(epiNodeLbl(nd))+(nd.aa?(' '+esc(nd.aa)):'')+'">'+esc(epiNodeLbl(nd))+'</th>';
      nodes.forEach(function(nd2,j){
        if(i===j){ h+='<td class="epimx-cell epimx-diag" title="'+esc(epiNodeLbl(nd))+' (self)"></td>'; return; }
        var v=cell(i,j);
        var t=esc(epiNodeLbl(nd))+' '+(v!=null&&v<0?'&#8646;':'&#8596;')+' '+esc(epiNodeLbl(nd2))+(v==null?': no shared series':': r = '+epiFmtR(v));
        h+='<td class="epimx-cell" style="background:'+epiColor(v)+'" title="'+t+'">'+((showVal&&v!=null)?('<span>'+epiFmtR(v)+'</span>'):'')+'</td>';
      });
      h+='</tr>';
    });
    h+='</tbody></table></div>';
    h+='<div class="epimx-scale"><span>discordant &#8722;1</span><i class="epimx-grad"></i><span>+1 concordant</span></div>';
    box.innerHTML=h;
  }
  function tsortVal(p,k){ if(k==='pair')return (p.geneA||'')+p.posA; if(k==='r')return p.r; if(k==='ar')return Math.abs(p.r); if(k==='n')return p.n; if(k==='p')return p.p; if(k==='q')return p.q; if(k==='tier')return ({strong:0,moderate:1,weak:2})[p.tier]; return p.direction; }
  function epiTable(list){
    var box=el('epi-table');
    var COLS=[['pair','Variant A &#8596; Variant B'],['direction','dynamics'],['ar','|r|'],['r','r'],['n','series'],['p','p'],['q','q (FDR)'],['tier','confidence']];
    var k=epiState.tsort.k, asc=epiState.tsort.asc;
    var rows=list.slice().sort(function(a,b){ var x=tsortVal(a,k),y=tsortVal(b,k),c; if(typeof x==='number'&&typeof y==='number')c=x-y; else c=String(x).localeCompare(String(y)); return asc?c:-c; });
    var shown=epiState.showAll?rows:rows.slice(0,12);
    var h='<div class="epitbl-top"><button class="dyn-btn" id="epidl" title="Download every reported pair as a TSV">'+icon('download')+'download pairs (TSV)</button><span class="dyn-count">'+rows.length+' pair(s)'+((!epiState.showAll&&rows.length>shown.length)?(' &#183; showing '+shown.length):'')+'</span></div>';
    h+='<div class="epitbl-wrap"><table class="epitbl"><thead><tr>'+COLS.map(function(c){return '<th data-k="'+c[0]+'">'+c[1]+(k===c[0]?(asc?icon('chevronUp','sort'):icon('chevronDown','sort')):'')+'</th>';}).join('')+'</tr></thead><tbody>';
    if(!rows.length){ h+='<tr><td colspan="'+COLS.length+'" class="c" style="padding:20px;text-align:center">no variant pair matches the current filter.</td></tr>'; }
    shown.forEach(function(p){
      var a='<b style="color:'+EPICOL.A+'">'+esc(p.geneA||'(intergenic)')+'</b>'+geneRvTag(p.geneA)+' '+String(p.posA).split(':').pop()+(p.aaA?(' '+aaDual(p.aaA,p.aaA_h37rv)):'');
      var b='<b style="color:'+EPICOL.B+'">'+esc(p.geneB||'(intergenic)')+'</b>'+geneRvTag(p.geneB)+' '+String(p.posB).split(':').pop()+(p.aaB?(' '+aaDual(p.aaB,p.aaB_h37rv)):'');
      h+='<tr><td>'+a+' <span class="epi-vs">'+(p.direction==='concordant'?'&#8596;':'&#8646;')+'</span> '+b+'</td>'+
        '<td>'+p.direction+'</td>'+
        '<td class="epitbl-r" style="color:'+(p.r>=0?'#2f8f5b':'#a24a8f')+'">'+(p.r>0?'+':'')+p.r.toFixed(2)+'</td>'+
        '<td>'+(p.r>0?'+':'')+p.r.toFixed(2)+'</td>'+
        '<td>'+p.n+(p.n>1?'&#8635;':'')+'</td>'+
        '<td>'+p.p.toFixed(3)+'</td>'+
        '<td>'+p.q.toFixed(3)+'</td>'+
        '<td><span class="epi-tier epi-'+p.tier+'">'+p.tier+'</span></td></tr>';
    });
    h+='</tbody></table></div>';
    box.innerHTML=h;
    Array.prototype.forEach.call(box.querySelectorAll('th[data-k]'),function(th){ th.onclick=function(){ var kk=th.getAttribute('data-k'); if(epiState.tsort.k===kk)epiState.tsort.asc=!epiState.tsort.asc; else{epiState.tsort.k=kk;epiState.tsort.asc=(kk==='pair'||kk==='direction'||kk==='tier');} epiTable(epiFilter()); }; });
    el('epidl').onclick=function(){
      var hdr=['geneA','posA','aa_A','geneB','posB','aa_B','dynamics','mean_r','n_series','r_min','r_max','perm_p','fdr_q','confidence'];
      var lines=[hdr.join('\t')];
      E.pairs.forEach(function(p){ lines.push([p.geneA,p.posA,p.aaA,p.geneB,p.posB,p.aaB,p.direction,p.r,p.n,p.rmin,p.rmax,p.p,p.q,p.tier].join('\t')); });
      dl(lines.join('\n')+'\n','epistasis_pairs.tsv','text/tab-separated-values');
    };
  }
  function draw(){
    var list=epiFilter();
    el('epicount').innerHTML=(epiState.view==='matrix')?((E.matrix?E.matrix.nodes.length:0)+' variant(s) in the matrix'):(list.length+' pair(s)'+((!epiState.showAll&&list.length>12)?(' &#183; showing 12'):'')+((epiState.dir==='all'&&epiState.conf==='all')?(' &#183; '+E.n_concordant+' concordant / '+E.n_discordant+' discordant &#183; '+E.n_strong+' strong'):''));
    el('epi-cards').style.display=epiState.view==='cards'?'':'none';
    el('epi-matrix').style.display=epiState.view==='matrix'?'':'none';
    el('epi-table').style.display=epiState.view==='table'?'':'none';
    if(epiState.view==='cards') epiCards(list);
    else if(epiState.view==='matrix') epiMatrix();
    else epiTable(list);
  }
  el('epiq').oninput=function(){ epiState.q=this.value; draw(); };
  Array.prototype.forEach.call(host.querySelectorAll('.epi-viewbtn'),function(b){ b.onclick=function(){ epiState.view=b.getAttribute('data-v'); Array.prototype.forEach.call(host.querySelectorAll('.epi-viewbtn'),function(x){x.className='epi-viewbtn'+(x.getAttribute('data-v')===epiState.view?' on':'');}); draw(); }; });
  Array.prototype.forEach.call(host.querySelectorAll('.epi-dirbtn'),function(b){ b.onclick=function(){ epiState.dir=b.getAttribute('data-d'); Array.prototype.forEach.call(host.querySelectorAll('.epi-dirbtn'),function(x){x.className='dyn-btn epi-dirbtn'+(x.getAttribute('data-d')===epiState.dir?' on':'');}); draw(); }; });
  Array.prototype.forEach.call(host.querySelectorAll('.epi-confbtn'),function(b){ b.onclick=function(){ epiState.conf=b.getAttribute('data-c'); Array.prototype.forEach.call(host.querySelectorAll('.epi-confbtn'),function(x){x.className='dyn-btn epi-confbtn'+(x.getAttribute('data-c')===epiState.conf?' on':'');}); draw(); }; });
  el('epir').oninput=function(){ epiState.minr=+this.value; el('epirv').textContent=epiState.minr.toFixed(2); draw(); };
  el('epizoom').oninput=function(){ epiZoom=+this.value; el('epi-cards').style.setProperty('--dyncw', epiZoom+'px'); };
  function epiSetAll(on){ epiState.showAll=on; var eb=el('epiShowAll'); if(eb){ eb.textContent=on?'show less':'show all'; eb.classList.toggle('on',on); } draw(); if(window.__syncSitesBtn)window.__syncSitesBtn(); }
  el('epiShowAll').onclick=function(){ epiSetAll(!epiState.showAll); };   // toggle: every reported pair <-> the top 12
  window.__epiSetAll=epiSetAll;
  if(epiState.showAll){ el('epiShowAll').textContent='show less'; el('epiShowAll').classList.add('on'); }
  draw();
}

