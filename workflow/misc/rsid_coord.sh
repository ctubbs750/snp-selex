IP_SNVS=$1
OP_RSID=$2

# awk '{print $11}' $IP_SNVS | while read rsid ;
# do
#   pos=$(curl -sX GET "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=$rsid&retmode=text&rettype=text" | sed 's/<\//\n/g' | grep -o -P '\<CHRPOS\>.{0,15}' | cut -f2 -d">") ;
# echo -e "${pos}" "${rsid}" | awk '{split($1,coord,":"); print "chr"coord[1],coord[2]-1, coord[2], $2}' | awk 'NF==4' >> $OP_RSID ;
# done ;

awk '{print $11}' $IP_SNVS | xargs -n 1 -P 24 -I {} bash -c '
export OP_RSID='"${OP_RSID}"'
rsid="{}"
pos=$(curl -sX GET "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=$rsid&retmode=text&rettype=text" | sed "s/<\//\n/g" | grep -o -P "\\<CHRPOS\\>.{0,15}" | cut -f2 -d">")
echo -e "${pos}" "${rsid}" | awk "{split(\$1,coord,\":\"); print \"chr\"coord[1],coord[2]-1, coord[2], \$2}" | awk "NF==4"  >> "${OP_RSID}"
'