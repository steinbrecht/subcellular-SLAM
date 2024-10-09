BEGIN { 
oldn=-1; 
oldm=-1; 
counter_exon=0;
counter_intron=0 } 

{ 
  if ((oldn!=$1) || (oldm!=$2)) { 
      if(oldn!=-1) {print oldn,oldm,counter_exon,counter_intron;} 
      oldn=$1; oldm=$2; counter_exon=0; counter_intron=0; 
  } 
  counter_exon=counter_exon+$3
  counter_intron=counter_intron+$4; 
}

#{print oldn, oldm, counter_exon, counter_intron}
