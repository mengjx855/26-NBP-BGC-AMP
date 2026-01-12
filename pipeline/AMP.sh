# tidy ref prot sequences
ls BGC_substring_pred/*pred | parallel -j 5 -k -q perl -e '$m=0;$n=$ARGV[0];while(<>){chomp;next if /ATCC 11775/; @s=split/\t/;for(1..$#s){$m+=$s[$_]};$m=$m/11;print "$n\t$_\t$m\n" if $m < 100; $m=0 }' {} > BGC.apex.result
perl -ne 'chomp; @s=split/\t/; $m=sprintf("%.6f",$s[-1]); $s[0]=~/\/(\S+?).faa/; print ">$1|$s[1]|$m\n$s[1]\n"' BGC.apex.result | seqkit rmdup -s -o BGC.apex.result.fa
seqkit fx2tab -n BGC.apex.result.fa | awk -F "|" '{print $0"\t"$3}' > BGC.apex.result.avgMIC
calcu_AMPs_similarity_self.py BGC.apex.result.fa BGC.apex.result.avgMIC BGC.apex.result.self
cut -f1 BGC.apex.result.self.unique | seqkit grep -f - BGC.apex.result.fa -o BGC.apex.result.self.unique.fa
calcu_AMPs_similarity_across.py BGC.apex.result.self.unique.fa amps.pep.ref.fa BGC.apex.result.across

# ref remove X [Unknown AA], B [D/N], Z [E/Q], J [I/L], U [硒半胱氨酸]
seqkit grep -s -v -r -p '[^ACDEFGHIKLMNPQRSTVWY*]' amps.pep.ref.fa > amps.pep.ref.removeXBZ.fa

# ref vs mine
seqkit grep -r -p 'DBAASPS.*' amps.pep.ref.fa | seqkit grep -s -v -r -p '[^ACDEFGHIKLMNPQRSTVWY*]' > DBAASPS.ref.fa
cat DBAASPS.ref.fa BGC.apex.result.self.unique.fa > all-to-all.2.fa

