cat Final_Blast_results.txt |while read line
	do  
	TAX=$(echo $line|awk '{print $14}') 
	NO_TAXONOMY='^[0-9]+$'
	 if [[ $TAX =~ $NO_TAXONOMY ]]; then 
		TAXID=$(echo $line|awk '{print $2}')
		TAXONOMY=$(/users/ybc502/edirect/efetch -db taxonomy -id $TAXID -format xml | xtract -pattern Taxon -first TaxId -element Taxon -block "*/Taxon" -unless Rank -equals "no rank" -tab " " -sep ":" -element Rank,ScientificName|perl -pe 's/(.+)(\s)order:(\w+) (.+)/$3/g')
			if [[ $TAXONOMY == "Lepidoptera" ]]; then
					 paste <(echo $line|awk 'NF{NF-=2}1') <(echo $TAXONOMY) <(echo $line|awk '{print $(NF-1)"\t"$NF}')
				else
					TAXONOMY=$(/users/ybc502/edirect/efetch -db taxonomy -id $TAXID -format xml | xtract -pattern Taxon -first TaxId -element Taxon -block "*/Taxon" -unless Rank -equals "no rank" -tab " " -sep ":" -element Rank,ScientificName||perl -pe 's/(.+)(\s)superkingdom:(\w+) (.+)/$3/g')					
					if [[ $TAXONOMY == "Bacteria" ]]; then
						paste <(echo $line|awk 'NF{NF-=2}1') <(echo $TAXONOMY) <(echo $line|awk '{print $(NF-1)"\t"$NF}')
					else
						TAXONOMY=$(/users/ybc502/edirect/efetch -db taxonomy -id $TAXID -format xml | xtract -pattern Taxon -first TaxId -element Taxon -block "*/Taxon" -unless Rank -equals "no rank" -tab " " -sep ":" -element Rank,ScientificName|perl -pe 's/(.+)(\s)class:(\w+) (.+)/$3/g')							
						paste <(echo $line|awk 'NF{NF-=2}1') <(echo $TAXONOMY) <(echo $line|awk '{print $(NF-1)"\t"$NF}')
					fi
				fi
			else
		echo $line
	fi
	done
