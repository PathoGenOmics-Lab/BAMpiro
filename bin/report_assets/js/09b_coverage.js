// ---- Deletions: the stretches no read covers, read against the rest of the cohort ----
// A stretch without reads is a deletion only when other samples on the same reference read it.
// When nobody does, it is a repeat or a part of the reference none of these genomes has, and
// those are listed apart rather than mixed in with what the samples lost.
var DEL_CLS={private:['private','only this sample lacks it; the other samples on its reference read it'],
  shared:['shared','several samples lack it; the others on the reference read it'],
  alone:['alone','the only sample on its reference read deeply enough, so nothing tells a deletion from a stretch the reference has and these genomes do not'],
  cohort:['nobody reads it','no sample reads it: a repeat, or a part of the reference none of these genomes has']};
var delView='del';
function delBins(r){   // a region's landscape bins: the landscape lays each sample's contigs end to end
  var C=R.coverage,nb=R.nbins||200,cs=(C.refs&&C.refs[r.ref])||[],off=0,tot=0,i;
  for(i=0;i<cs.length;i++){if(cs[i][0]==r.contig)off=tot;tot+=cs[i][1];}
  if(!tot)tot=R.genome_len||1;
  function b(p){return Math.max(0,Math.min(nb-1,Math.floor((off+p)*nb/tot)));}
  return [b(r.start),b(r.end)];}
function delWhere(x){   // 'contig:start-end', with the reference first when a run has several and it is not the contig
  var multi=Object.keys((R.coverage&&R.coverage.refs)||{}).length>1;
  return (multi&&x.ref!=x.contig?esc(x.ref)+' &#183; ':'')+esc(x.contig)+':'+x.start.toLocaleString('en-US')+'&#8211;'+x.end.toLocaleString('en-US');}
function delMatch(r,q){
  if(!q)return true;
  if((r.contig+':'+r.start+'-'+r.end).toLowerCase().indexOf(q)>=0)return true;
  if(r.genes.some(function(g){return g.toLowerCase().indexOf(q)>=0;}))return true;
  return r.samples.some(function(x){return String(x[0]).toLowerCase().indexOf(q)>=0;});}
function renderDeletions(){
  var host=el('del_body'); if(!host)return;
  var C=R.coverage;
  if(!(C&&C.regions)){host.innerHTML='';return;}
  var q=(st.delq||'').toLowerCase(), all=C.regions;
  var cnt={private:0,shared:0,alone:0,cohort:0}; all.forEach(function(r){cnt[r.cls]=(cnt[r.cls]||0)+1;});
  var nDel=cnt.private+cnt.shared+cnt.alone;
  var regs=all.filter(function(r){return ((delView=='cohort')==(r.cls=='cohort'))&&delMatch(r,q);});
  var big=el('deletionsPanel')&&el('deletionsPanel').classList.contains('expanded'), cap=q?400:(big?300:40);
  var lede='<b>'+nDel+'</b> '+_plural(nDel,'stretch','stretches')+' of '+C.min_len+' bp or more that some samples do not read and others do: '+
    '<b>'+cnt.private+'</b> private to one sample, <b>'+cnt.shared+'</b> shared'+(cnt.alone?(', '+cnt.alone+' in a sample alone on its reference'):'')+'. '+
    '<b>'+cnt.cohort+'</b> '+_plural(cnt.cohort,'stretch','stretches')+' no sample reads '+_plural(cnt.cohort,'is','are')+' set apart: a repeat or a part of the reference these genomes lack is nobody&#39;s deletion. '+
    'Measured on the '+C.n_assessed+' samples read at a median depth of '+C.min_depth+'&#215; or more'+
    (C.not_assessed&&C.not_assessed.length?('; '+C.not_assessed.length+' thinner '+_plural(C.not_assessed.length,'sample was','samples were')+' left out, since stretches without reads turn up there by chance'):'')+'.';
  var seg='<span class="seg" id="delseg"><button data-v="del"'+(delView=='del'?' class="on"':'')+'>Deletions '+nDel+'</button>'+
    '<button data-v="cohort"'+(delView=='cohort'?' class="on"':'')+'>Nobody reads '+cnt.cohort+'</button></span>';
  var rows=regs.slice(0,cap).map(function(r){
    var bb=delBins(r), names=r.samples.map(function(x){return x[0];});
    var who=_list(names.map(function(s){return esc(s);}),4);
    var tipS=r.samples.map(function(x){return x[0]+': '+x[1].toLocaleString('en-US')+'-'+x[2].toLocaleString('en-US');}).join('\n');
    var genes=r.genes.length?_list(r.genes.map(function(g){return esc(g)+geneRvTag(g);}),6):'<span class="na">none</span>';
    return '<tr data-b0="'+bb[0]+'" data-b1="'+bb[1]+'" data-id="'+esc(r.id)+'">'+
      '<td class="s" style="text-align:left">'+delWhere(r)+'</td>'+
      '<td>'+r.len.toLocaleString('en-US')+'</td>'+
      '<td title="'+esc(tipS)+'">'+r.n+' of '+r.of+'</td>'+
      '<td style="text-align:left" title="'+esc(tipS)+'">'+who+'</td>'+
      '<td style="text-align:left">'+genes+'</td>'+
      '<td style="text-align:left" title="'+esc(DEL_CLS[r.cls][1])+'">'+esc(DEL_CLS[r.cls][0])+'</td></tr>';}).join('');
  var head='<tr><th class="s" style="text-align:left">Stretch</th><th>bp</th><th>Samples</th><th style="text-align:left">Which</th><th style="text-align:left">Genes it removes</th><th style="text-align:left">Class</th></tr>';
  var empty='<tr><td colspan="6" class="na" style="text-align:left;padding:8px">'+(q?('nothing matches &quot;'+esc(q)+'&quot;.'):(delView=='cohort'?'every stretch without reads is read by some sample.':'no sample lacks a stretch the others read.'))+'</td></tr>';
  var more=regs.length>cap?('<div class="hot-note">Showing '+cap+' of '+regs.length+'. Search, open the panel full, or download the table for the rest.</div>'):'';
  var gl=C.genes_lost||[], gtab='';
  if(gl.length&&delView=='del'){
    gtab='<div class="dsub" style="margin:14px 16px 4px">Genes some samples do not read <span style="font-weight:400;color:#94a3b8">(under half of the gene read, where most samples on the reference read nearly all of it)</span></div>'+
      '<div class="gtable" style="max-height:32vh"><table><thead><tr><th class="s" style="text-align:left">Gene</th><th style="text-align:left">Position</th><th>Samples</th><th style="text-align:left">Which (share of the gene read)</th></tr></thead><tbody>'+
      gl.filter(function(g){return !q||g.gene.toLowerCase().indexOf(q)>=0||g.samples.some(function(x){return String(x[0]).toLowerCase().indexOf(q)>=0;});}).slice(0,cap).map(function(g){
        return '<tr><td class="s" style="text-align:left">'+esc(g.gene)+geneRvTag(g.gene)+'</td><td style="text-align:left">'+delWhere(g)+'</td>'+
          '<td>'+g.samples.length+' of '+g.of+'</td>'+
          '<td style="text-align:left">'+_list(g.samples.map(function(x){return esc(x[0])+' ('+Math.round(x[1]*100)+'%)';}),6)+'</td></tr>';}).join('')+
      '</tbody></table></div>';}
  host.innerHTML='<div class="hot-note" style="margin:10px 16px 6px">'+lede+'</div>'+
    '<div style="display:flex;gap:10px;align-items:center;margin:0 16px 8px;flex-wrap:wrap">'+seg+
    '<button class="dyn-btn" id="deldl" title="Download every stretch as a TSV, with each sample&#39;s own coordinates">'+icon('download')+'download (TSV)</button></div>'+
    '<div class="gtable" style="max-height:'+(big?'70vh':'44vh')+'"><table id="deltable"><thead>'+head+'</thead><tbody>'+(rows||empty)+'</tbody></table></div>'+more+gtab;
  Array.prototype.forEach.call(host.querySelectorAll('#delseg button'),function(b){b.onclick=function(){delView=b.getAttribute('data-v');renderDeletions();};});
  var dlb=el('deldl'); if(dlb)dlb.onclick=function(){
    var lines=['region\treference\tcontig\tstart\tend\tlength\tclass\tn_samples\tn_assessed\tsamples\tgenes'];
    all.forEach(function(r){lines.push([r.id,r.ref,r.contig,r.start,r.end,r.len,r.cls,r.n,r.of,
      r.samples.map(function(x){return x[0]+':'+x[1]+'-'+x[2];}).join(','),r.genes.join(',')].join('\t'));});
    dl(lines.join('\n')+'\n','deletions.tsv','text/tab-separated-values');};
  Array.prototype.forEach.call(host.querySelectorAll('#deltable tbody tr[data-b0]'),function(tr){tr.onclick=function(){
    var b0=+tr.getAttribute('data-b0'),b1=+tr.getAttribute('data-b1'),nb=R.nbins||200,pad=Math.max(3,Math.round((b1-b0)*0.6)+2);
    if(gtrackHas('del'))st.gtrack='del';
    st.geneMark={b0:b0,b1:b1,name:tr.getAttribute('data-id')||''}; st.gzoom={b0:Math.max(0,b0-pad),b1:Math.min(nb-1,b1+pad)};
    Array.prototype.forEach.call(document.querySelectorAll('#gtrack button'),function(x){x.classList.toggle('on',x.getAttribute('data-gt')==st.gtrack);});
    renderGenome(); el('genome').scrollIntoView();};});
}
