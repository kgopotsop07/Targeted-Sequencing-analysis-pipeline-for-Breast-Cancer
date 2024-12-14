#!/bin/bash

usage(){
	printf "============help===============\n"
	printf "\e[4mget_analysis_databases.sh:\n\e[0m"
	printf "	get_analysis_databases.sh -r <reference type>\n"
	printf "\e[4m\nPipeline option:\n\e[0m"
	printf "	-r:.................. specify the reference database type to be downloaded <hg19 or hg38>\n"
	printf "get_analysis_databases.sh -h | --help\n\n" && exit
	}

optstring=":r:h:"
while getopts ${optstring} arg; do
	case ${arg} in
		r) r=${OPTARG} ; ((r == "hg19" || r == "hg38")) ;;
		h) h=${OPTARG} usage;;
		?) ;;
		*) ;;
	esac
done

if [ -z $1 ] || [ $1 = "--help" ] || [ $1 = "-h" ]; then
	usage && exit
fi

## ----------- input variables ----------------
DB_TYPE=${r}
variant_reference_files="${PWD}/reference_files/variant_calling_reference_files/$DB_TYPE"
reference_file_path="${PWD}/reference_files/human_genome_reference_file/$DB_TYPE"
annotation_reference_files="${PWD}/reference_files/annotationt_reference_files/$DB_TYPE"
echo "$variant_reference_files"
## ----------- links to variant calling databases ------------------
## ---- Genome Reference Consortium 37 (hg19)
if [ "$DB_TYPE" == "hg19" ]; then
KNOWN_INDELS="https://console.cloud.google.com/storage/browser/_details/gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.known_indels_20120518.vcf?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&inv=1&invt=Abj7fQ"
KNOWN_SNPS="https://console.cloud.google.com/storage/browser/_details/gcp-public-data--broad-references/hg19/v0/1000G_phase1.snps.high_confidence.b37.vcf.gz?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&inv=1&invt=Abj7fQ"
KNOWN_MILLS_GOLD_STANDARD_INDELS="https://console.cloud.google.com/storage/browser/_details/gcp-public-data--broad-references/hg19/v0/Mills_and_1000G_gold_standard.indels.b37.sites.vcf?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&inv=1&invt=Abj7fQ"
PON="gs://gatk-best-practices/somatic-b37/Mutect2-exome-panel.vcf"
DBSNP="https://console.cloud.google.com/storage/browser/_details/gcp-public-data--broad-references/hg19/v0/dbsnp_135.b37.vcf.gz?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))&inv=1&invt=Abj7fQ"

elif [ "$DB_TYPE" == "hg38" ]; then
## ---- Genome Reference Consortium 38 (hg38) ----
KNOWN_INDELS="https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))"
KNOWN_SNPS="https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz?pageState=(%22StorageObjectListTable%22:(%22f%22:%22%255B%255D%22))"
KNOWN_MILLS_GOLD_STANDARD_INDELS="https://console.cloud.google.com/storage/browser/_details/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PON="gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz"

fi

get_GCh19(){
	echo " => downloading gatk resource bundle known indels for ${DB_TYPE}"
	 wget -c $KNOWN_INDELS -o '${variant_reference_files}/hg19_indels.knownsites.vcf.gz'
	 wget -c $KNOWN_SNPS -o '${variant_reference_files}/hg19_snps.knownsites.vcf.gz'
	 wget -c $KNOWN_MILLS_GOLD_STANDARD_INDELS -o '${hg19_variant_reference_files}/hg19_Mills_indels.knownsites.vcf'
	 wget -c $PON -o '${reference_file_path}/hg19_somatic1000g_pon.b37.vcf'
	 wget -c DBSNP -o '${reference_file_path}/hg19_af.only.gnomad.vcf'

	echo " => fetching annotation database for annovar $DB_TYPE"
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar clinvar_20240611 annotation_reference_files
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar refGene annotation_reference_files
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar avsnp147 annotation_reference_files
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar cosmic70 annotation_reference_files

}

get_GCh38(){
	echo " => downloading gatk resource bundle known indels for ${DB_TYPE}"
	 wget -c $KNOWN_INDELS -o '${variant_reference_files}/hg38_indels.knownsites.vcf.gz'
	 wget -c $KNOWN_SNPS -o '$variant_reference_files/hg38_snps.knownsites.vcf.gz'
	 wget -c $KNOWN_MILLS_GOLD_STANDARD_INDELS -o '$hg38_variant_reference_files/hg38_Mills_indels.knownsites.vcf'
	 wget -c DBSNP -o '$reference_file_path/hg38_af.only.gnomad.vcf'

	echo " => fetching annotation database for annovar $DB_TYPE"
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar clinvar_20240611 annotation_reference_files
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar refGene annotation_reference_files
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar avsnp147 annotation_reference_files
	 annotate_variation.pl -buildver $DB_TYPE -downdb -webfrom annovar cosmic70 annotation_reference_files
}

if [ "$DB_TYPE" == "hg19" ]; then
	get_GCh19
elif [ "$DB_TYPE" == "hg38" ]; then
	get_GCh38
else
	usage && exit
fi
