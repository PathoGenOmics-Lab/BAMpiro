function snpAfColor(af){return 'rgba(31,120,180,'+(0.16+af*0.8).toFixed(2)+')';}
var SNPMX_PAL=['#bcd0ea','#f3d1b0','#c3e0c9','#f0c4cf','#d6c9ec','#b8e0dd','#eadfb0','#dfe4ea','#f2c4c4','#cdd1a8','#e6c3e0','#b9d6ee'];
// dark-theme categorical palette for the SNP-matrix metadata rows (same hues, muted/dark so they don't glare as light bands)
var SNPMX_PAL_DARK=['#2f4a6b','#6b4c2f','#2f5a40','#6b3a4a','#4a3a6b','#2f5a55','#5c4a2f','#3a4450','#6b3838','#4c4c2f','#5a2f5a','#35506b'];
var snpmxFilter={};
var snpmxAll=false;   // false = show the top MAXR most-shared sites; true = virtualized scroll over every site
function renderSnpMatrix(){
  var host=el('snpmx_body'), sec=el('snpmatrix'); if(!host)return;
  var M=R.snp_matrix;
  if(!(M&&M.rows&&M.rows.length)){ if(sec)sec.style.display='none'; var nv=el('nav-snpmx'); if(nv)nv.style.display='none'; return; }
  if(sec)sec.style.display='';
  var samples=M.samples, MAXR=400;
  var meta=(R.sample_meta&&R.sample_meta.fields&&R.sample_meta.fields.length)?R.sample_meta:null;
  var metaMaps={}, metaVals={};
  var PAL=isDark()?SNPMX_PAL_DARK:SNPMX_PAL;
  if(meta){ meta.fields.forEach(function(f){ var m={},vals=[],k=0; samples.forEach(function(s){var v=(meta.rows[s]||{})[f]; if(v&&!(v in m)){m[v]=PAL[k%PAL.length];k++;vals.push(v);}}); metaMaps[f]=m; metaVals[f]=vals.sort(); }); }
  // a field is "the lineage field" only if every (non-NA) value is a known lineage -> reuse the canonical
  // report palette (mycolorsTB / linColor) for it; all other fields get the neutral pastels. Keying on the
  // field (not just the value) stops an unrelated column whose value happens to equal a lineage label from
  // being mis-coloured.
  var linField={};
  if(meta){ meta.fields.forEach(function(f){ var vs=(metaVals[f]||[]).filter(function(v){return v&&v!=='NA'&&v!=='.'&&v!=='-';}); if(vs.length&&vs.every(function(v){return LINCOL[v];})) linField[f]=1; }); }
  function metaColor(f,v){ if(v&&linField[f]&&LINCOL[v]) return LINCOL[v]; return (v&&metaMaps[f]&&metaMaps[f][v])?metaMaps[f][v]:TH.cellnull; }
  function metaText(bg){
    if(bg&&bg.charAt(0)==='#'){ var h=bg.length===4?('#'+bg.charAt(1)+bg.charAt(1)+bg.charAt(2)+bg.charAt(2)+bg.charAt(3)+bg.charAt(3)):bg;
      var L=(0.299*parseInt(h.substr(1,2),16)+0.587*parseInt(h.substr(3,2),16)+0.114*parseInt(h.substr(5,2),16))/255; return L<0.62?'#fff':'#1c2b3a'; }
    var m=/hsl\(\s*[\d.]+\s*,\s*[\d.]+%\s*,\s*([\d.]+)%/i.exec(bg||'');   // algorithmic lineage fallback is hsl(h,58%,52%)
    if(m) return (+m[1])<62?'#fff':'#1c2b3a';
    return '#1c2b3a';
  }
  function visIdx(){ var idx=[]; samples.forEach(function(s,i){ var ok=true; if(meta){ for(var f in snpmxFilter){ if(snpmxFilter[f] && (meta.rows[s]||{})[f]!==snpmxFilter[f]){ ok=false; break; } } } if(ok) idx.push(i); }); return idx; }
  var filterUI=meta?('<div class="snpmx-filters"><span class="snpmx-flabel">filter columns:</span>'+
    meta.fields.map(function(f){ return '<label class="snpmx-fsel">'+esc(f)+' <select data-f="'+esc(f)+'"><option value="">all</option>'+
      metaVals[f].map(function(v){return '<option value="'+esc(v)+'"'+(snpmxFilter[f]===v?' selected':'')+'>'+esc(v)+'</option>';}).join('')+'</select></label>'; }).join('')+
    '<button class="dyn-btn" id="snpmxfclear">clear</button></div>'):'';
  host.innerHTML=
    '<div class="snpmx-controls">'+
      '<input id="snpmxq" class="dyn-search" type="search" placeholder="filter by gene / position / amino acid...">'+
      '<label class="snpmx-toggle"><input type="checkbox" id="snpmxdp" checked> show depth</label>'+
      '<button class="dyn-btn" id="snpmxdl" title="Download the full matrix (all samples) as a wide TSV">'+icon('download')+'download matrix (TSV)</button>'+
      '<button class="dyn-btn snpmx-allbtn" id="snpmxallbtn" style="display:none" title="Toggle between the top most-shared sites and a scrollable view of every site"></button>'+
      '<span class="dyn-count" id="snpmxcount"></span></div>'+
    filterUI+
    (meta?('<div class="snpmx-metanote">column levels from the samplesheet: '+meta.fields.map(function(f){return '<b>'+esc(f)+'</b>';}).join(' &#183; ')+' &#183; hover a header cell for its value</div>'):'')+
    '<div class="snpmx-wrap" id="snpmxwrap"><table class="snpmx" id="snpmxtable"></table></div>';
  function draw(){
    var q=(el('snpmxq').value||'').toLowerCase(), showDP=el('snpmxdp').checked;
    var vi=visIdx();
    var rows=M.rows.filter(function(r){
      if(q && !((r.gene&&r.gene.toLowerCase().indexOf(q)>=0)||String(r.pos).indexOf(q)>=0||(r.aa&&r.aa.toLowerCase().indexOf(q)>=0))) return false;
      return vi.some(function(i){return r.cells[i];});   // only SNPs seen in the visible columns
    });
    var nfilt=0; for(var kf in snpmxFilter){ if(snpmxFilter[kf]) nfilt++; }
    var total=rows.length, capped=total>MAXR, ncol=1+vi.length;
    el('snpmxcount').innerHTML=total+' SNP site(s) &#215; '+vi.length+' sample(s)'+(nfilt?' (filtered)':'')+
      (capped?(' &#183; showing <b>'+(snpmxAll?('all '+total):(MAXR+' of '+total))+'</b>'):'')+(M.truncated?' &#183; full matrix in the TSV':'');
    var allbtn=el('snpmxallbtn');   // a real button toggles between the top sites and the full virtualized scroll
    if(allbtn){ if(capped){ allbtn.style.display=''; allbtn.classList.toggle('on',snpmxAll);
        allbtn.innerHTML=snpmxAll?('show top '+MAXR):('show all '+total+' &#8595;');
        allbtn.onclick=function(){ snpmxAll=!snpmxAll; el('snpmxwrap').scrollTop=0; draw(); if(window.__syncSitesBtn)window.__syncSitesBtn(); }; }
      else { allbtn.style.display='none'; } }
    var mh=22, nf=meta?meta.fields.length:0;
    var metaRows=meta?meta.fields.map(function(f,k){
      return '<tr>'+'<th class="snpmx-info snpmx-metalabel" style="top:'+(k*mh)+'px">'+esc(f)+'</th>'+
        vi.map(function(i){var s=samples[i],v=(meta.rows[s]||{})[f]||'',bg=metaColor(f,v); return '<th class="snpmx-metacell" style="top:'+(k*mh)+'px;background:'+bg+';color:'+metaText(bg)+'" title="'+esc(f)+': '+esc(v||'-')+'">'+esc(v)+'</th>';}).join('')+'</tr>';
    }).join(''):'';
    var stop=nf*mh;
    var nameRow='<tr><th class="snpmx-info snpmx-corner" style="top:'+stop+'px">SNP '+esc(M.reference?('('+M.reference+')'):'')+'</th>'+
      vi.map(function(i){var s=samples[i];return '<th class="snpmx-hcell" style="top:'+stop+'px" title="'+esc(s)+'"><span class="snpmx-h">'+esc(s)+'</span></th>';}).join('')+'</tr>';
    function rowHTML(r){
      var lbl='<b>'+esc(r.gene||r.contig)+'</b>'+geneRvTag(r.gene)+' '+r.pos+' '+esc(r.ref)+'&#8594;'+esc(r.alt)+(r.aa?(' <span class="snpmx-aa">'+aaDual(r.aa,r.aa_h37rv)+'</span>'):'');
      var cells=vi.map(function(i){var c=r.cells[i], s=samples[i];
        if(!c) return '<td class="snpmx-cell snpmx-empty" title="'+esc(s)+' - not called"></td>';
        var afTxt=c[0].toFixed(2).replace(/^0/,'').replace(/^1\.00$/,'1');
        var dpTxt=(showDP&&c[1]!=null)?('<span class="snpmx-dp">'+c[1]+'</span>'):'';
        return '<td class="snpmx-cell'+(showDP?' wdp':'')+'" style="background:'+snpAfColor(c[0])+'" title="'+esc(s)+'  AF='+c[0].toFixed(3)+(c[1]!=null?('  DP='+c[1]):'')+'"><span class="snpmx-af">'+afTxt+'</span>'+dpTxt+'</td>';
      }).join('');
      return '<tr><td class="snpmx-info">'+lbl+'</td>'+cells+'</tr>';
    }
    var tbl=el('snpmxtable'), thead='<thead>'+metaRows+nameRow+'</thead>', wrap=el('snpmxwrap');
    wrap.onscroll=null;
    if(snpmxAll&&capped){
      // virtualized scroll: render only the rows in view, with top/bottom spacer rows keeping the scrollbar honest
      tbl.innerHTML=thead+'<tbody id="snpmxbody"></tbody>';
      var tb=el('snpmxbody'), RH=27, lastStart=-1;
      var win=function(force){
        var stp=wrap.scrollTop, vh=wrap.clientHeight||500;
        var per=Math.ceil(vh/RH)+10, start=Math.max(0,Math.floor(stp/RH)-5), end=Math.min(total,start+per);
        if(!force&&start===lastStart)return; lastStart=start;
        var topH=start*RH, botH=(total-end)*RH, h='';
        if(topH>0)h+='<tr class="snpmx-spacer"><td colspan="'+ncol+'" style="height:'+topH+'px"></td></tr>';
        for(var i=start;i<end;i++)h+=rowHTML(rows[i]);
        if(botH>0)h+='<tr class="snpmx-spacer"><td colspan="'+ncol+'" style="height:'+botH+'px"></td></tr>';
        tb.innerHTML=h;
      };
      win(true);
      var probe=tb.querySelector('tr:not(.snpmx-spacer)'); if(probe){ var ph=probe.offsetHeight; if(ph>6&&Math.abs(ph-RH)>1){ RH=ph; lastStart=-1; win(true); } }
      // repaint the window straight from the scroll event; win() cheaply no-ops until the start row actually changes
      wrap.onscroll=function(){ win(false); };
    } else {
      tbl.innerHTML=thead+'<tbody>'+rows.slice(0,MAXR).map(rowHTML).join('')+'</tbody>';
    }
  }
  el('snpmxq').oninput=function(){clearTimeout(_mxdb);_mxdb=setTimeout(draw,160);};
  el('snpmxdp').onchange=draw;
  if(meta){
    Array.prototype.forEach.call(host.querySelectorAll('.snpmx-filters select'),function(sel){ sel.onchange=function(){ var f=sel.getAttribute('data-f'); if(sel.value) snpmxFilter[f]=sel.value; else delete snpmxFilter[f]; draw(); }; });
    el('snpmxfclear').onclick=function(){ snpmxFilter={}; Array.prototype.forEach.call(host.querySelectorAll('.snpmx-filters select'),function(s){s.value='';}); draw(); };
  }
  el('snpmxdl').onclick=function(){
    var hdr=['reference','contig','pos','ref_allele','alt_allele','gene','effect','aa_change'];
    samples.forEach(function(s){hdr.push(s+'|AF');hdr.push(s+'|DP');});
    var lines=[hdr.join('\t')];
    if(meta){ meta.fields.forEach(function(f){ var row=['# '+f,'','','','','','','']; samples.forEach(function(s){row.push((meta.rows[s]||{})[f]||'');row.push('');}); lines.push(row.join('\t')); }); }
    M.rows.forEach(function(r){var row=[M.reference||r.contig,r.contig,r.pos,r.ref,r.alt,r.gene,r.eff,r.aa];
      samples.forEach(function(s,i){var c=r.cells[i]; if(c){row.push(c[0].toFixed(4));row.push(c[1]==null?'':c[1]);}else{row.push('');row.push('');}});
      lines.push(row.join('\t'));});
    dl(lines.join('\n')+'\n','snp_matrix.tsv','text/tab-separated-values');
  };
  window.__snpmxDraw=draw;   // let the header 'all sites' shortcut repaint the matrix without a full re-render
  draw();
  if(window.__syncSitesBtn)window.__syncSitesBtn();
}

